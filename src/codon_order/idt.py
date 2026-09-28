"""IDT SciTools Plus API client — eBlock synthesis-complexity screening.

This is the *manufacturability gate* for synthesis prep, not a codon optimizer:
codon optimization runs locally in ``codon_optimize.py`` (DnaChisel), and IDT is
asked only "would you actually make this?".  The number returned is the same
**Total Complexity Score** the ordering site shows, so the manual
"re-optimise the hard ones" loop becomes a pipeline step.

Endpoint (confirmed from the Baker lab's SAPP/DMX release, ``JB/idt.py``)::

    POST https://www.idtdna.com/Restapi/v1/Complexities/ScreenEBlockSequences
    body: [{"Name": "...", "Sequence": "ATG..."}, ...]

The response is one list per submitted sequence: an **empty** list means the
sequence raised no issue at all (score 0.0); otherwise each entry carries a
``Score`` and the sum over entries is the total.  Lower is better; the project
gate is ``< 7`` (see ``docs/decisions/0001-codon-optimization.md``).

This screens *manufacturability only*.  IDT's own order-time biosecurity
screening is separate and unaffected by anything here.

Credentials are read from a git-ignored ``.env`` (never hard-coded, never
logged); see ``.env.example``. Two credential sets are needed:
  * the API client you created (``IDT_CLIENT_ID`` / ``IDT_CLIENT_SECRET``), and
  * your IDT account login (``IDT_USERNAME`` / ``IDT_PASSWORD``), used at token
    time (OAuth2 *password* grant).

Usage::

    python idt_complexity.py --check-auth                  # validate credentials
    python idt_complexity.py --sequence ATGAAA...          # score one sequence
    codon-order-complexity --csv batch1.csv
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import requests

TOKEN_URL = "https://www.idtdna.com/Identityserver/connect/token"
COMPLEXITY_URL = "https://www.idtdna.com/Restapi/v1/Complexities/ScreenEBlockSequences"

# IDT caps how many sequences one screening call accepts; keep batches modest.
MAX_BATCH = 20

REQUIRED_VARS = ("IDT_CLIENT_ID", "IDT_CLIENT_SECRET", "IDT_USERNAME", "IDT_PASSWORD")


def load_env(env_path: Path = Path(".env")) -> dict[str, str]:
    """Parse a simple KEY=VALUE ``.env`` into a dict (no external dependency).

    Values already present in the real environment win, so you can also export
    the vars instead of using the file.
    """
    vals: dict[str, str] = {}
    if env_path.exists():
        for line in env_path.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            vals[k.strip()] = v.strip()
    # real environment overrides file
    for k in REQUIRED_VARS:
        if os.environ.get(k):
            vals[k] = os.environ[k]
    missing = [k for k in REQUIRED_VARS if not vals.get(k)]
    if missing:
        raise SystemExit(f"missing credentials in .env / environment: {', '.join(missing)}")
    return vals


def credentials_available(env_path: Path = Path(".env")) -> bool:
    """True if all four credentials are resolvable, without raising."""
    try:
        load_env(env_path)
    except SystemExit:
        return False
    return True


def get_token(env: dict[str, str], session: requests.Session) -> dict:
    """Exchange credentials for a short-lived bearer token (never logged)."""
    resp = session.post(
        TOKEN_URL,
        auth=(env["IDT_CLIENT_ID"], env["IDT_CLIENT_SECRET"]),  # HTTP Basic
        data={
            "grant_type": "password",
            "username": env["IDT_USERNAME"],
            "password": env["IDT_PASSWORD"],
            "scope": "test",
        },
        timeout=30,
    )
    if resp.status_code != 200:
        # IDT error bodies describe the failure but do not echo the password.
        raise SystemExit(
            f"token request failed: HTTP {resp.status_code}\n{resp.text[:500]}"
        )
    tok = resp.json()
    if not tok.get("access_token"):
        raise SystemExit(f"no access_token in response: {list(tok)}")
    return tok


class ComplexityClient:
    """Bearer-token session that scores eBlock sequences on demand.

    The token is fetched lazily on first use and refreshed shortly before it
    expires, so a long optimisation run (many sequences x several retries) needs
    no special handling by the caller.
    """

    def __init__(self, env_path: Path = Path(".env"), refresh_margin: int = 60) -> None:
        self._env = load_env(env_path)
        self._session = requests.Session()
        self._margin = refresh_margin
        self._token: str | None = None
        self._expires_at = 0.0
        self.calls = 0

    def _bearer(self) -> str:
        if self._token is None or time.time() > self._expires_at - self._margin:
            tok = get_token(self._env, self._session)
            self._token = tok["access_token"]
            self._expires_at = time.time() + float(tok.get("expires_in", 3600))
        return self._token

    def screen(self, named_sequences: list[tuple[str, str]]) -> list[float]:
        """Return one total complexity score per (name, sequence) pair."""
        scores: list[float] = []
        for start in range(0, len(named_sequences), MAX_BATCH):
            batch = named_sequences[start:start + MAX_BATCH]
            payload = [{"Name": n, "Sequence": s} for n, s in batch]
            resp = self._session.post(
                COMPLEXITY_URL,
                headers={"Content-Type": "application/json",
                         "Authorization": f"Bearer {self._bearer()}"},
                data=json.dumps(payload),
                timeout=120,
            )
            self.calls += 1
            if resp.status_code != 200:
                raise RuntimeError(
                    f"complexity screening failed: HTTP {resp.status_code}\n"
                    f"{resp.text[:500]}"
                )
            issues = resp.json()
            if len(issues) != len(batch):
                raise RuntimeError(
                    f"response length {len(issues)} != {len(batch)} submitted"
                )
            # An empty issue list means "nothing flagged" -> score 0.
            scores.extend(
                float(sum(entry.get("Score", 0.0) for entry in per_seq or []))
                for per_seq in issues
            )
        return scores

    def score(self, sequence: str, name: str = "eBlock") -> float:
        """Total complexity score for a single sequence."""
        return self.screen([(name, sequence)])[0]

    def details(self, sequence: str, name: str = "eBlock") -> list[dict]:
        """Raw per-issue records for one sequence (for explaining a failure)."""
        resp = self._session.post(
            COMPLEXITY_URL,
            headers={"Content-Type": "application/json",
                     "Authorization": f"Bearer {self._bearer()}"},
            data=json.dumps([{"Name": name, "Sequence": sequence}]),
            timeout=120,
        )
        self.calls += 1
        if resp.status_code != 200:
            raise RuntimeError(
                f"complexity screening failed: HTTP {resp.status_code}\n{resp.text[:500]}"
            )
        return resp.json()[0] or []


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check-auth", action="store_true",
                    help="Only validate credentials by fetching a token, then exit.")
    ap.add_argument("--sequence", help="Score a single DNA sequence and print its issues.")
    ap.add_argument("--csv", help="Score every row of a codon_optimized_*.csv "
                                  "(expects 'name' and 'sequence' columns).")
    ap.add_argument("--env", default=".env", help="Path to the credentials file.")
    args = ap.parse_args()

    if args.check_auth:
        env = load_env(Path(args.env))
        with requests.Session() as session:
            tok = get_token(env, session)
        print(f"✅ token acquired (expires_in={tok.get('expires_in')}s, "
              f"token_type={tok.get('token_type')}, "
              f"length={len(tok['access_token'])})")
        return

    client = ComplexityClient(Path(args.env))

    if args.sequence:
        seq = args.sequence.strip().upper()
        issues = client.details(seq)
        total = sum(i.get("Score", 0.0) for i in issues)
        print(f"{len(seq)} bp -> total complexity score {total:.1f}")
        for i in issues:
            print(f"  {i.get('Score', 0.0):>5.1f}  {i.get('Name', '?')}: "
                  f"{i.get('DisplayText', '')}")
        return

    if args.csv:
        import pandas as pd

        df = pd.read_csv(args.csv)
        pairs = list(zip(df["name"].astype(str), df["sequence"].astype(str)))
        for (name, _), score in zip(pairs, client.screen(pairs)):
            print(f"  {score:>6.1f}  {name}")
        return

    ap.error("nothing to do: pass --check-auth, --sequence or --csv")


if __name__ == "__main__":
    main()

"""CLI karpos_surfex — orchestration de la chaîne SURFEX OFFLINE.

    python -m karpos_surfex run --workdir run/ --ecoclimap MY_RUN/ECOCLIMAP \
                                --steps pgd,prep,offline
    python -m karpos_surfex where          # localise les exécutables
"""

from __future__ import annotations

import argparse
import sys

from . import __version__, driver, orchestrate
from ._surfex import SurfexError


def main(argv=None):
    p = argparse.ArgumentParser(prog="karpos_surfex", description=__doc__)
    p.add_argument("--version", action="version", version=f"karpos-surfex {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="enchaîne PGD→PREP→OFFLINE dans un répertoire")
    r.add_argument("--workdir", required=True, help="répertoire de run (namelists + forcing)")
    r.add_argument("--ecoclimap", default=None, help="répertoire des covers ECOCLIMAP (*.bin)")
    r.add_argument("--steps", default="pgd,prep,offline", help="étapes CSV (pgd,prep,offline,soda)")

    w = sub.add_parser("where", help="localise les exécutables SURFEX")
    w.add_argument("--prog", default=None, help="un maître précis (PGD/PREP/OFFLINE/SODA)")

    args = p.parse_args(argv)

    try:
        if args.cmd == "run":
            steps = tuple(s.strip() for s in args.steps.split(",") if s.strip())
            res = orchestrate.run_chain(args.workdir, steps=steps, ecoclimap_dir=args.ecoclimap)
            print(f"OK — étapes {res['steps']} dans {res['workdir']}")
            if res["ecoclimap"]:
                print(f"     ECOCLIMAP liés : {res['ecoclimap']}")
            print(f"     sorties : {res['outputs']}")
        elif args.cmd == "where":
            progs = [args.prog.upper()] if args.prog else list(orchestrate.driver.MASTERS)
            for prog in progs:
                print(f"{prog:8s} → {driver.find_exe(prog)}")
    except SurfexError as e:
        print(f"ERREUR : {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

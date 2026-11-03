"""Command-line entry point: verify, recompute, or render saved data."""
from pathlib import Path
import argparse,json

def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    sub.add_parser('verify',help='Derive the exact tensor identities; no output files')
    r=sub.add_parser('recompute',help='Recompute formulas and figure data in a NEW directory');r.add_argument('--output',type=Path,required=True)
    f=sub.add_parser('render',help='Render saved data; no symbolic computation');f.add_argument('--data',type=Path,default=Path(__file__).resolve().parents[1]/'data/paper03ec_v4');f.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='verify':
        from paper_results.derive import derive
        from paper_results.extensions import derive_extensions
        results=[derive(),derive_extensions()]
        print(json.dumps({'exact_identities':sum(len(r['identity_residuals']) for r in results),'all_zero':all(v=='0' for r in results for v in r['identity_residuals'].values())}))
    elif args.command=='recompute':
        from paper_results.data import generate
        generate(args.output);print(args.output)
    else:
        from paper_results.figures import render
        for path in render(args.data,args.output):print(path)

if __name__=='__main__':main()

import argparse
from pathlib import Path
import sys

from .config import load_config, read_genes
from .http import HttpClient, RecordedSession
from .pipeline import annotate_stage, associate_stage, fetch_stage, report_stage
from .storage import Store


def main(argv=None):
    parser = argparse.ArgumentParser(description="Gene annotation lab with public API data and synthetic sample data")
    parser.add_argument("command", choices=["fetch", "annotate", "associate", "report", "all"])
    parser.add_argument("--config", type=Path, default=Path("config.yaml"))
    parser.add_argument("--genes", type=Path, default=Path("data/genes.txt"))
    parser.add_argument("--offline", action="store_true", help="Use recorded public HTTP responses without network access")
    args = parser.parse_args(argv)
    store, client = None, None
    try:
        config = load_config(args.config)
        symbols = read_genes(args.genes)
        store = Store(config["output"])
        client = HttpClient(config, RecordedSession(config["recordings"]) if args.offline else None)
        if args.command in {"fetch", "all"}:
            fetch_stage(symbols, client, config, store)
        if args.command in {"annotate", "all"}:
            annotate_stage(client, config, store)
        if args.command in {"associate", "all"}:
            associate_stage(config, store)
        if args.command in {"report", "all"}:
            print(report_stage(symbols, config, store))
        print(f"Completed {args.command}. Outputs: {config['output']}")
        return 0
    except (ValueError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    finally:
        if client is not None:
            client.close()
        if store is not None:
            store.close()


if __name__ == "__main__":
    raise SystemExit(main())

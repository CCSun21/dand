import sys
from dandelion import dandelion_prep, dandelion_sample, dandelion_refine

def main():
    if len(sys.argv) < 2:
        print("Usage: dand [prep|sample|refine] [options]")
        sys.exit(1)

    command = sys.argv[1]
    # Remove the 'dand' and the subcommand from sys.argv
    sys.argv = [sys.argv[0]] + sys.argv[2:]

    if command == "prep":
        dandelion_prep.main()
    elif command == "sample":
        dandelion_sample.main()
    elif command == "refine":
        dandelion_refine.main()
    else:
        print(f"Unknown command: {command}")
        print("Available commands: prep, sample, refine")
        sys.exit(1)

if __name__ == "__main__":
    main()
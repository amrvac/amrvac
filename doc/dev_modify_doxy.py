from pathlib import Path


def main():
    doxyfile = Path(__file__).parent.resolve() / "Doxyfile"
    assert doxyfile.is_file()
    with open(doxyfile, "r") as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith("PROJECT_NUMBER"):
            lines[i] = 'PROJECT_NUMBER = "3.1"\n'
        if line.startswith("PROJECT_BRIEF"):
            brief = lines[i].split("=")[1].replace('"', "").strip()
            lines[i] = f'PROJECT_BRIEF = "{brief} (development version)"\n'
    with open(doxyfile, "w") as f:
        f.writelines(lines)


if __name__ == "__main__":
    main()

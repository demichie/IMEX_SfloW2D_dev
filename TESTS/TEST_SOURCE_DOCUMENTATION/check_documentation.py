"""Check Doxygen coverage of production Fortran modules and procedures.

Run with Python 3; no compiler, Doxygen installation or third-party packages
are required. This is a signature/header consistency check, not a replacement
for Doxygen rendering or a scientific review of the descriptions.
"""
from pathlib import Path
import re

START = re.compile(
    r'^\s*(?:(?:pure|elemental|recursive|impure)\s+'
    r'|(?:real|complex|integer|logical|character|double precision)'
    r'(?:\s*\([^)]*\))?\s+)*(subroutine|function)\s+(\w+)'
    r'\s*(?:\(([^)]*)\))?', re.I)
MODULE = re.compile(r'^\s*(module|program)\s+(\w+)\s*$', re.I)
END = re.compile(r'^\s*end\s*(subroutine|function|module|program)\b', re.I)

def code_part(line):
    """Remove a Fortran comment without removing exclamation marks in strings."""
    quote = None
    i = 0
    while i < len(line):
        c = line[i]
        if quote:
            if c == quote:
                if i + 1 < len(line) and line[i + 1] == quote:
                    i += 1
                else:
                    quote = None
        elif c in "\"'":
            quote = c
        elif c == '!':
            return line[:i]
        i += 1
    return line

def split_list(text):
    """Split comma-separated declarations without splitting array dimensions."""
    items, start, depth = [], 0, 0
    for i, c in enumerate(text):
        if c == '(':
            depth += 1
        elif c == ')':
            depth -= 1
        elif c == ',' and depth == 0:
            items.append(text[start:i].strip())
            start = i + 1
    items.append(text[start:].strip())
    return [i for i in items if i]

def statements(lines):
    """Join continued free-form statements and retain their source locations."""
    pending, start = '', 0
    for i, line in enumerate(lines):
        code = code_part(line).strip()
        if not code:
            continue
        if not pending:
            start = i
        pending += ' ' + code.lstrip('&').rstrip('&').strip()
        if code.endswith('&'):
            continue
        yield start, i, pending.strip()
        pending = ''
    assert not pending, 'Unterminated continuation'

def entities(path):
    """Collect procedure/module headers and declared dummy-argument directions."""
    lines = path.read_text().splitlines()
    found, stack = [], []
    for first, last, statement in statements(lines):
        match = START.match(statement) or MODULE.match(statement)
        if match:
            kind, name = match.group(1).lower(), match.group(2)
            args = []
            if kind in ('subroutine', 'function') and match.group(3):
                args = split_list(match.group(3))
            begin = first
            while begin and (not lines[begin-1].strip()
                             or lines[begin-1].lstrip().startswith('!')):
                begin -= 1
            entity = dict(kind=kind, name=name, args=args, start=first,
                          header_start=begin, end=None, intents={},
                          declarations={}, header=lines[begin:first])
            found.append(entity)
            stack.append(entity)
        elif END.match(statement):
            kind = END.match(statement).group(1).lower()
            assert stack and stack[-1]['kind'] == kind, (path, first, stack)
            stack.pop()['end'] = last
        elif stack and '::' in statement:
            attributes, variables = statement.split('::', 1)
            intent = re.search(r'\bintent\s*\(\s*(inout|out|in)\s*\)',
                               attributes, re.I)
            dummy_names = {arg.lower() for arg in stack[-1]['args']}
            for variable in split_list(variables):
                name_match = re.match(r'(\w+)', variable)
                if name_match and name_match.group(1).lower() in dummy_names:
                    name = name_match.group(1).lower()
                    stack[-1]['declarations'][name] = statement
                    if intent:
                        stack[-1]['intents'][name] = intent.group(1).lower()
    assert not stack, (path, stack)
    return found


# These expressions deliberately cover the project's free-form declarations,
# including continued signatures, internal procedures and interface prototypes.
PARAM = re.compile(
    r"[\\@]param[ \t]*\[(in(?:[ \t]*,?[ \t]*out)?|out)\]"
    r"[ \t]+(\w+)[ \t]+(\S.*)", re.I)
BRIEF = re.compile(r"[\\@]brief[ \t]+\S", re.I)
RETURN = re.compile(r"[\\@]returns?[ \t]+\S", re.I)


def validate_header(entity):
    """Return signature/header discrepancies, without modifying the source."""
    header = "\n".join(entity["header"])
    errors = []
    if not BRIEF.search(header):
        errors.append("missing Doxygen brief")
    params = {}
    for direction, name, _description in PARAM.findall(header):
        name = name.lower()
        if name in params:
            errors.append(f"duplicate parameter documentation: {name}")
        params[name] = re.sub(r"[\s,]", "", direction.lower())
    expected = {arg.lower() for arg in entity["args"]}
    for name in sorted(expected - params.keys()):
        errors.append(f"undocumented parameter: {name}")
    for name in sorted(params.keys() - expected):
        errors.append(f"documented parameter absent from signature: {name}")
    for name, intent in entity["intents"].items():
        if name in params and params[name] != intent:
            errors.append(f"{name}: documented {params[name]}, declared INTENT({intent})")
    if entity["kind"] == "function" and not RETURN.search(header):
        errors.append("missing function return documentation")
    return errors


def executable_lines(source):
    """Preserve executable code, preprocessing and compiler/OpenMP directives."""
    result = []
    for line in source.splitlines():
        stripped = line.strip()
        if stripped.startswith(("!$", "!DEC$", "!DIR$", "!GCC$")):
            result.append(stripped)
        else:
            code = code_part(line).strip()
            if code:
                result.append(code)
    return result


def main():
    import argparse
    from collections import Counter
    import subprocess

    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument(
        "--source-dir", type=Path,
        default=Path(__file__).resolve().parents[2] / "src",
        help="production Fortran source directory (default: repository src)")
    cli.add_argument(
        "--check-against-git", metavar="REF",
        help="also require executable lines and directives to match a Git revision")
    args = cli.parse_args()
    source_dir = args.source_dir.resolve()
    paths = sorted(source_dir.glob("*.f90"))
    if not paths:
        cli.error(f"no .f90 sources found in {source_dir}")

    counts = Counter()
    errors = []
    for path in paths:
        for entity in entities(path):
            counts[entity["kind"]] += 1
            for message in validate_header(entity):
                errors.append(f"{path.name}:{entity['start']+1}: {entity['name']}: {message}")

    fragments = sorted(source_dir.glob("*.inc"))
    for path in fragments:
        header = path.read_text()
        if not re.search(r"[\\@]file\s+" + re.escape(path.name), header):
            errors.append(f"{path.name}: missing Doxygen file block")
        if not BRIEF.search(header):
            errors.append(f"{path.name}: missing Doxygen brief")

    if args.check_against_git:
        repo = source_dir.parent
        for path in paths + fragments:
            reference = subprocess.check_output(
                ["git", "-C", str(repo), "show",
                 f"{args.check_against_git}:{path.relative_to(repo).as_posix()}"],
                text=True)
            if executable_lines(reference) != executable_lines(path.read_text()):
                errors.append(f"{path.name}: executable code/directives differ from {args.check_against_git}")

    if errors:
        print("\n".join(errors))
        return 1
    summary = ", ".join(f"{count} {kind}{'s' if count != 1 else ''}"
                        for kind, count in sorted(counts.items()))
    print(f"PASS: {len(paths)} source files, {summary}, {len(fragments)} include fragments")
    if args.check_against_git:
        print(f"PASS: executable code and directives unchanged from {args.check_against_git}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

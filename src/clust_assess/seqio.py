
from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, Iterator

QUERY_PREFIX = "Q__"
_BAD = re.compile(r"[^A-Za-z0-9._-]+")


def sanitise(raw: str) -> str:
    s = _BAD.sub("_", raw.strip()).strip("_.")
    return s or "unnamed"


def query_id(raw: str) -> str:
    return QUERY_PREFIX + sanitise(raw)


def is_query(pid: str) -> bool:
    return pid.startswith(QUERY_PREFIX)


def input_name(pid: str) -> str:
    return pid[len(QUERY_PREFIX):] if is_query(pid) else pid


def read_fasta(path: str | Path) -> Iterator[tuple[str, str, str]]:
    name, desc, chunks = None, "", []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, desc, "".join(chunks)
                head = line[1:].strip()
                name, _, desc = head.partition(" ")
                chunks = []
            elif line.strip():
                chunks.append(line.strip())
    if name is not None:
        yield name, desc, "".join(chunks)


def write_fasta(path: str | Path, records: Iterable[tuple[str, str]], wrap: int = 80) -> None:
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), wrap):
                fh.write(seq[i:i + wrap] + "\n")


def header(pid: str, circular: bool) -> str:
    return f"{pid} circular=true" if circular else pid

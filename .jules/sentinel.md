## 2025-05-14 - Buffer Overflows and Memory Leaks in lacepr
**Vulnerability:** Multiple buffer overflows via `sscanf`, `strcpy`, and `strcat`, plus a major memory leak in the BAM processing loop.
**Learning:** The C component `lacepr` used unsafe string functions without length checks when parsing recalibration files and BAM tags. It also allocated memory for strings on every read/record without freeing them.
**Prevention:** Always use width limits in `sscanf` (e.g., `%255s`), use `strncpy` with explicit null-termination, and avoid redundant heap allocations in high-frequency loops by pointing to existing buffers or reusing memory.

## 2025-05-15 - Command Injection via two-argument open() in Perl
**Vulnerability:** Use of two-argument `open()` with user-supplied filenames in `lacer.pl`.
**Learning:** In Perl, the two-argument `open(FH, $filename)` can execute commands if `$filename` starts with a pipe (`|`). This is a classic command injection vector.
**Prevention:** Always use the three-argument form `open(my $fh, '<', $filename)` and lexical filehandles. This explicitly separates the open mode from the filename, preventing injection.

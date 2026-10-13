import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import make_depend as cd


class TestLogicalLines(unittest.TestCase):
    def test_fixed_form_continuation_and_comments(self):
        # fixed-form basics: a column-1 'C' is a whole-line comment,
        # and a non-blank column 6 continues the previous statement,
        # so the two source lines become one logical line -- reported
        # at the line it started on, which is what warnings cite.
        src = (
            "C this is a comment, ignored\n"
            "      call foo(a,\n"
            "     &          b, c)\n"
            "      x = 1\n"
        )
        p = Path(self.id() + ".F")
        p.write_text(src)
        self.addCleanup(p.unlink)
        lines = list(cd.iter_logical_lines(p))
        self.assertEqual(len(lines), 2)
        self.assertIn("call foo(a,", lines[0].text)
        self.assertIn("b, c)", lines[0].text)
        self.assertEqual(lines[0].first_lineno, 2)

    def test_fixed_form_cpp_directive_is_not_column_stripped(self):
        # regression: a #include in a .F file is not a Fortran
        # statement, so column 1 is the '#' itself, not a label
        # field -- the fixed-form column-6 stripping must not apply
        # to it (caught during Phase 4, since column 6 of "#include"
        # happens to look like a continuation marker).
        src = '#include "missing.h"\n      x=1\n'
        p = Path(self.id() + ".F")
        p.write_text(src)
        self.addCleanup(p.unlink)
        lines = list(cd.iter_logical_lines(p))
        self.assertEqual(len(lines), 2)
        self.assertEqual(lines[0].text, '#include "missing.h"')
        self.assertEqual(lines[1].text, "x=1")

    def test_fixed_form_inline_bang_comment_is_stripped(self):
        # regression: real .F files in this tree (main.F, step.F) use
        # trailing "! comment" after a statement, not just the
        # classic column-1 whole-line convention. Found via Phase 7's
        # differential run against the real tree: "USE xios ! XIOS
        # module" was being misparsed as UnresolvedUse("xios")
        # because the comment was never stripped.
        src = "      USE xios           ! XIOS module\n"
        p = Path(self.id() + ".F")
        p.write_text(src)
        self.addCleanup(p.unlink)
        lines = list(cd.iter_logical_lines(p))
        self.assertEqual(len(lines), 1)
        self.assertNotIn("XIOS module", lines[0].text)

    def test_fixed_form_tab_indented_statement(self):
        # regression: the GNU/DEC "tab form" extension, which
        # cross_matrix.F handled explicitly (it tested for char(9) in
        # the label field) but this tool did not. A leading TAB ends
        # the label field, so the statement starts right after it --
        # slicing at column 6 instead eats the first five characters
        # and loses the dependency entirely. Real occurrences in this
        # tree: online_interpolate_bulk.F:319, init_scalars.F:2406.
        src = "\tuse mymod\n\tcall foo(a)\n"
        p = Path(self.id() + ".F")
        p.write_text(src)
        self.addCleanup(p.unlink)
        lines = [ll.text for ll in cd.iter_logical_lines(p)]
        self.assertEqual(lines, ["use mymod", "call foo(a)"])

    def test_fixed_form_tab_then_digit_is_a_continuation(self):
        # in tab form a nonzero digit right after the TAB is the
        # continuation marker, not part of the statement
        src = "      x = 1 +\n\t12\n"
        p = Path(self.id() + ".F")
        p.write_text(src)
        self.addCleanup(p.unlink)
        lines = [ll.text for ll in cd.iter_logical_lines(p)]
        self.assertEqual(lines, ["x = 1 +2"])

    def test_fixed_form_tab_in_column_six_starts_a_statement(self):
        # checked against gfortran: a TAB in column 6 ends the label
        # field and starts a new statement, it is NOT read as the
        # column-6 continuation marker
        src = "      x = 1 +\n     \t2\n"
        p = Path(self.id() + ".F")
        p.write_text(src)
        self.addCleanup(p.unlink)
        lines = [ll.text for ll in cd.iter_logical_lines(p)]
        self.assertEqual(lines, ["x = 1 +", "2"])

    def test_label_before_tab_is_not_part_of_the_statement(self):
        # tab form again, this time with a statement label in front of
        # the TAB: the label belongs to the label field, not to the
        # statement
        src = "10\tcontinue\n"
        p = Path(self.id() + ".F")
        p.write_text(src)
        self.addCleanup(p.unlink)
        lines = [ll.text for ll in cd.iter_logical_lines(p)]
        self.assertEqual(lines, ["continue"])

    def test_c_style_block_comment_is_stripped(self):
        # regression: this tree's headers use C-style /* ... */
        # comments, sometimes spanning multiple lines, e.g.
        # cppdefs_dev.h's "/* use options from NBQ ... */" -- found
        # via Phase 7's differential run: prose inside such a block
        # was being misread as a real USE statement, since nothing
        # stripped it. (_strip_block_comment blanks the comment text
        # rather than deleting it, so that fixed-form column positions
        # survive; this fixture is free-form and does not exercise
        # that part.)
        src = (
            "      x = 1\n"
            "/*\n"
            "   use options from NBQ Open boundary conditions\n"
            "*/\n"
            "      y = 2\n"
        )
        p = Path(self.id() + ".F90")
        p.write_text(src)
        self.addCleanup(p.unlink)
        lines = [ll.text for ll in cd.iter_logical_lines(p)]
        joined = " ".join(lines)
        self.assertNotIn("use options", joined)
        self.assertIn("x = 1", joined)
        self.assertIn("y = 2", joined)

    def test_block_comment_inside_single_line(self):
        # a /* ... */ that opens and closes on one line must not
        # swallow the rest of the statement with it
        src = "      x = /* inline block */ 1\n"
        p = Path(self.id() + ".F90")
        p.write_text(src)
        self.addCleanup(p.unlink)
        lines = list(cd.iter_logical_lines(p))
        self.assertEqual(len(lines), 1)
        self.assertNotIn("inline block", lines[0].text)
        self.assertIn("x =", lines[0].text)
        self.assertIn("1", lines[0].text)

    def test_free_form_continuation_and_bang_comment(self):
        # the free-form half of iter_logical_lines' suffix-based
        # dispatch: '!' is the only comment marker, and a trailing '&'
        # continues onto the next line regardless of any column
        src = (
            "  x = 1  ! inline comment\n"
            "  call foo(a, &\n"
            "           b)\n"
        )
        p = Path(self.id() + ".F90")
        p.write_text(src)
        self.addCleanup(p.unlink)
        lines = list(cd.iter_logical_lines(p))
        self.assertEqual(len(lines), 2)
        self.assertNotIn("inline comment", lines[0].text)
        self.assertIn("call foo(a,", lines[1].text)
        self.assertIn("b)", lines[1].text)


class TestTokenize(unittest.TestCase):
    def test_identifiers_and_eol(self):
        # every stream ends with exactly one EOL token; the USE
        # grammar leans on it to tell a bare "use modb" from a
        # truncated one, without special-casing the end of the list
        toks = cd.tokenize("use modb")
        self.assertEqual(
            [t.kind for t in toks],
            [cd.TokKind.IDENT, cd.TokKind.IDENT, cd.TokKind.EOL],
        )
        self.assertEqual([t.value for t in toks[:2]], ["use", "modb"])

    def test_string_literal_single_and_double_quoted(self):
        # both quote styles have to produce a STRING carrying the
        # unquoted filename: Fortran INCLUDE uses '...', CPP #include
        # uses "...", and this tree contains both
        toks = cd.tokenize("include 'mpif.h'")
        self.assertEqual(toks[1].kind, cd.TokKind.STRING)
        self.assertEqual(toks[1].value, "mpif.h")
        toks2 = cd.tokenize('#include "param.h"')
        strtok = [t for t in toks2 if t.kind is cd.TokKind.STRING][0]
        self.assertEqual(strtok.value, "param.h")

    def test_doubled_quote_is_an_escaped_quote(self):
        # '' inside a quoted string is one literal quote, not the end
        # of the string. No filename in this tree needs it; this is
        # lexer correctness, so the scanner cannot lose its place.
        toks = cd.tokenize("include 'it''s.h'")
        self.assertEqual(toks[1].value, "it's.h")

    def test_comma_coloncolon_and_hash(self):
        # the punctuation the USE grammar has to distinguish -- ','
        # (attribute or only-clause), '::' and a leading '#' -- gets
        # its own token kind rather than falling into OTHER
        toks = cd.tokenize("use, intrinsic :: iso_c_binding")
        kinds = {t.kind for t in toks}
        self.assertIn(cd.TokKind.COMMA, kinds)
        self.assertIn(cd.TokKind.COLONCOLON, kinds)
        self.assertEqual(cd.tokenize("#include")[0].kind, cd.TokKind.HASH)


class TestLooksRelevant(unittest.TestCase):
    """The pre-filter that lets scan_file skip tokenizing the ~96% of
    logical lines that cannot be a dependency-relevant statement."""

    def test_accepts_keywords_and_hash(self):
        # every statement kind the parser can accept has to survive
        # the filter, in any case and with any leading whitespace
        for text in ("use modb", "  Use modb", "include 'x.h'",
                     "module foo", '#include "x.h"'):
            self.assertTrue(cd.looks_relevant(text), text)

    def test_rejects_everything_else(self):
        # ...and the filter must key off the whole first word, not a
        # prefix: "used" and "user_type" both begin with "use" and a
        # naive substring check would wrongly let them through
        for text in ("x = 1", "call foo(a, b)", "used = .true.",
                     "user_type % field = 1", ""):
            self.assertFalse(cd.looks_relevant(text), text)


class TestParseStatement(unittest.TestCase):
    def _parse(self, text):
        return cd.parse_statement(cd.tokenize(text))

    def test_mixed_case_use_is_recognized(self):
        # Fortran is case-insensitive, but cross_matrix.F matched the
        # literal 'USE '/'use ' only, so the real "Use Agrif_Grids" in
        # Agrif2Model.F and zoom.F was invisible to it. Module names
        # are normalised to lower case for the registry lookup.
        stmt = self._parse("Use Agrif_Grids")
        self.assertIsInstance(stmt, cd.UseStmt)
        self.assertEqual(stmt.module, "agrif_grids")

    def test_use_comma_no_space(self):
        # the F2003 attribute form with every optional space removed.
        # Absent from this tree today, but it is standard Fortran and
        # cross_matrix.F could not parse it at any spacing, since it
        # required a literal 'use ' with a trailing blank.
        stmt = self._parse("use,intrinsic::iso_c_binding")
        self.assertIsInstance(stmt, cd.UseStmt)
        self.assertEqual(stmt.module, "iso_c_binding")

    def test_use_only_clause_keeps_module_name(self):
        # the dependency is the module, not the only-list. This is the
        # dominant shape in CROCO ("use croco_namelist, only: dt"),
        # and the one cross_matrix.F needed a strip-a-trailing-comma
        # hack to half-handle.
        stmt = self._parse("use modb, only: foo")
        self.assertIsInstance(stmt, cd.UseStmt)
        self.assertEqual(stmt.module, "modb")

    def test_bare_fortran_include(self):
        # Fortran INCLUDE rather than CPP #include. Real in this tree
        # -- MPI_Setup.F, MessPass2D.F and MessPass3D.F all do
        # "include 'mpif.h'" -- and never recognised by
        # cross_matrix.F, which only ever looked for #include "...".
        stmt = self._parse("INCLUDE 'mpif.h'")
        self.assertIsInstance(stmt, cd.IncludeStmt)
        self.assertEqual(stmt.filename, "mpif.h")

    def test_cpp_include(self):
        # the CPP form, which is how CROCO pulls in every one of its
        # headers and the main thing Make.depend has to track
        stmt = self._parse('#include "param.h"')
        self.assertIsInstance(stmt, cd.IncludeStmt)
        self.assertEqual(stmt.filename, "param.h")

    def test_module_procedure_is_not_a_definition(self):
        # "module procedure" inside an interface block is not a module
        # definition; accepting it would register a module literally
        # named "procedure" (real: exchange_device_host.F:18 and :24)
        self.assertIsNone(self._parse("MODULE PROCEDURE foo"))

    def test_module_definition(self):
        # definitions are what populate the registry, so that a USE
        # resolves to the file that really defines the module instead
        # of to a guessed filename. Stored lower-cased, to match the
        # normalisation done on the USE side.
        stmt = self._parse("module ModB")
        self.assertIsInstance(stmt, cd.ModuleDefStmt)
        self.assertEqual(stmt.name, "modb")

    def test_unresolved_use_is_flagged_not_dropped(self):
        # a USE target that is not a plain identifier (a macro-built
        # module name, say) is reported rather than quietly ignored,
        # so it can be checked by hand. Not observed in this tree --
        # if this warning ever fires for real, go and look.
        stmt = self._parse("use MOD_NAME_MACRO(x)")
        self.assertIsInstance(stmt, cd.UnresolvedUse)
        self.assertEqual(stmt.raw, "MOD_NAME_MACRO")


class TestFilterNeverStricterThanParser(unittest.TestCase):
    """The pre-filter is an optimisation, so it must not change the
    result. One-directional on purpose: the filter may pass a line the
    parser then rejects, but never the reverse."""

    CASES = [
        "use modb", "      Use Agrif_Grids",
        "use,intrinsic::iso_c_binding", "use modb, only: foo",
        "INCLUDE 'mpif.h'", '#include "param.h"',
        "module ModB", "MODULE PROCEDURE foo",
        "use MOD_NAME_MACRO(x)", "#define FOO 1", "#ifdef BAR",
        "x = 1", "call foo(a, b)", "used = .true.",
        "user_type % field = 1", 'call foo("include me")',
        "implicit none", "",
    ]

    def test_filter_is_never_stricter_than_the_parser(self):
        # checked over every statement shape the parser knows about,
        # plus near-misses chosen to tempt a sloppy filter
        for text in self.CASES:
            stmt = cd.parse_statement(cd.tokenize(text))
            if stmt is not None:
                self.assertTrue(
                    cd.looks_relevant(text),
                    f"filter would drop a real statement: {text!r} -> {stmt}",
                )


class TestCatalogAndRegistry(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(self.id())
        self.tmp.mkdir(exist_ok=True)
        self.addCleanup(lambda: [p.unlink() for p in self.tmp.glob("*")] and self.tmp.rmdir())

    def _write(self, name, text):
        p = self.tmp / name
        p.write_text(text)
        return p

    def test_missing_include_warns_and_is_excluded(self):
        # a header that is not on disk is warned about and catalogued
        # as found=False, rather than aborting the scan: AGRIF
        # generates some of its headers at build time, and
        # netcdf.inc/mpif.h live outside the tree altogether
        a = self._write("a.F", '#include "missing.h"\n      x=1\n')
        catalog, warnings = cd.build_catalog([a])
        self.assertTrue(any("missing.h" in w for w in warnings))
        self.assertIn("missing.h", catalog)

    def test_registry_resolves_module_case_insensitively(self):
        # "module ModB" and "use MODB" name the same module. Matching
        # case-sensitively is what made cross_matrix.F invent four
        # separate files for one AGRIF module (AGRIF_Util.F,
        # Agrif_Util.F, AGRIF_UTIL.F, Agrif_UTIL.F), none of which
        # exist.
        self._write("modb.F90", "      module ModB\n      end module\n")
        a = self._write("a.F90", "      use MODB\n")
        catalog, _ = cd.build_catalog([a, self.tmp / "modb.F90"])
        registry, _ = cd.build_module_registry(catalog)
        self.assertEqual(registry["modb"], "modb.F90")

    def test_duplicate_module_definition_keeps_first_and_warns(self):
        # two files claiming one module name is ambiguous, so resolve
        # it deterministically (first wins) and say so rather than
        # letting dict order decide. No occurrence in the tree today;
        # this pins the behaviour for when one appears.
        self._write("m1.F90", "      module dup\n      end module\n")
        self._write("m2.F90", "      module dup\n      end module\n")
        catalog, _ = cd.build_catalog([self.tmp / "m1.F90", self.tmp / "m2.F90"])
        registry, warnings = cd.build_module_registry(catalog)
        self.assertEqual(registry["dup"], "m1.F90")
        self.assertTrue(any("dup" in w for w in warnings))

    def test_fast_filter_toggle_agrees_with_full_scan(self):
        # raised 2026-09-25: keep the filter removable, not just
        # trusted -- and prove the two paths agree on a fixture
        # that exercises every statement kind at once.
        a = self._write(
            "a.F90",
            "      use modb\n"
            "      Use modb\n"
            "      include 'x.h'\n"
            "      x = 1\n",
        )
        self._write("modb.F90", "      module modb\n      end module\n")
        self._write("x.h", "      integer :: n\n")
        fast_catalog, fast_warn = cd.build_catalog([a], fast_filter=True)
        full_catalog, full_warn = cd.build_catalog([a], fast_filter=False)
        self.assertEqual(fast_catalog["a.F90"].uses, full_catalog["a.F90"].uses)
        self.assertEqual(fast_catalog["a.F90"].includes, full_catalog["a.F90"].includes)
        self.assertEqual(fast_warn, full_warn)


class TestResolveDependencies(unittest.TestCase):
    def test_missing_include_is_excluded_from_resolved_deps(self):
        # regression, found via a real `make` run on the BASIN test
        # case: a sequential (non-MPI) build has no mpif.h, correctly
        # warned-and-excluded by build_catalog -- but it was still
        # being emitted as a literal prerequisite here, so `make`
        # failed outright with "No rule to make target 'mpif.h'".
        # cross_matrix.F's own EXCLUDE actually drops such files from
        # the emitted dependency list, not just warns about them.
        catalog = {
            "a.F": cd.FileInfo(Path("a.F"), {"missing.h"}, set(), set(), True),
            "missing.h": cd.FileInfo(
                Path("missing.h"), set(), set(), set(), False, found=False
            ),
        }
        deps = cd.resolve_dependencies(catalog, {})
        self.assertEqual(deps["a.F"], set())

    def test_header_closure_is_transitive(self):
        # headers are not build targets, so nothing else in the
        # makefile will ever pull in a.F's indirect headers: the full
        # #include closure has to be listed on a.F's own line here
        catalog = {
            "a.F": cd.FileInfo(Path("a.F"), {"b.h"}, set(), set(), True),
            "b.h": cd.FileInfo(Path("b.h"), {"c.h"}, set(), set(), False),
            "c.h": cd.FileInfo(Path("c.h"), set(), set(), set(), False),
        }
        deps = cd.resolve_dependencies(catalog, {})
        self.assertEqual(deps["a.F"], {"b.h", "c.h"})

    def test_mutual_header_inclusion_terminates(self):
        # headers in this tree do include each other back. The closure
        # has to terminate, and must not end up making a.F its own
        # prerequisite, which make would read as a circular dependency.
        catalog = {
            "a.F": cd.FileInfo(Path("a.F"), {"b.h"}, set(), set(), True),
            "b.h": cd.FileInfo(Path("b.h"), {"a.F"}, set(), set(), False),
        }
        deps = cd.resolve_dependencies(catalog, {})
        self.assertEqual(deps["a.F"], {"b.h"})

    def test_use_inside_included_header_bubbles_up(self):
        # a USE sitting inside an #include'd header is a USE by the
        # including file, since the header is spliced in before the
        # compiler sees it -- so a.F, not b.h, depends on modc
        catalog = {
            "a.F": cd.FileInfo(Path("a.F"), {"b.h"}, set(), set(), True),
            "b.h": cd.FileInfo(Path("b.h"), set(), {"modc"}, set(), False),
            "modc.F90": cd.FileInfo(Path("modc.F90"), set(), set(), {"modc"}, False),
        }
        registry = {"modc": "modc.F90"}
        deps = cd.resolve_dependencies(catalog, registry)
        # module dependencies are converted to their .o form here,
        # in resolve_dependencies, not later in the writer -- see
        # test_included_raw_source_file_keeps_its_own_name for why
        self.assertEqual(deps["a.F"], {"b.h", "modc.o"})

    def test_included_raw_source_file_keeps_its_own_name(self):
        # regression, found via Phase 7's differential run against
        # the real tree: CROCO commonly #includes a raw .F file as
        # spliced-in source rather than USEing a module (e.g.
        # biology.F #include-ing bio_BioEBUS.F). Such a file is not
        # a separate build target, so it must never be rewritten to
        # .o the way a real USE'd module is -- both can end in .F,
        # so the distinction must be made by *how* the dependency
        # was discovered (#include vs USE), not by suffix alone.
        catalog = {
            "biology.F": cd.FileInfo(
                Path("biology.F"), {"bio_BioEBUS.F"}, set(), set(), True
            ),
            "bio_BioEBUS.F": cd.FileInfo(
                Path("bio_BioEBUS.F"), set(), set(), set(), False
            ),
        }
        deps = cd.resolve_dependencies(catalog, {})
        self.assertEqual(deps["biology.F"], {"bio_BioEBUS.F"})


class TestFindUnresolvedUses(unittest.TestCase):
    def test_module_defined_nowhere_is_reported(self):
        # raised in review: silence here is how renaming the file that
        # defines a module drops its dependency edge without a trace.
        # A NOTE rather than a WARNING, because the usual cause is a
        # genuinely external module, as xios is here.
        catalog = {
            "a.F": cd.FileInfo(Path("a.F"), set(), {"xios"}, set(), True),
        }
        notes = cd.find_unresolved_uses(catalog, {})
        self.assertEqual(len(notes), 1)
        self.assertIn("xios", notes[0])
        self.assertIn("a.F", notes[0])
        self.assertTrue(notes[0].startswith("NOTE:"))

    def test_resolved_module_is_not_reported(self):
        # the ordinary case has to stay quiet, or the notes are noise
        # and nobody reads the one that matters
        catalog = {
            "a.F": cd.FileInfo(Path("a.F"), set(), {"modc"}, set(), True),
            "modc.F90": cd.FileInfo(Path("modc.F90"), set(), set(), {"modc"}, False),
        }
        self.assertEqual(cd.find_unresolved_uses(catalog, {"modc": "modc.F90"}), [])

    def test_use_inside_an_included_header_is_attributed_to_the_includer(self):
        # same bubbling-up rule as resolve_dependencies, so the note
        # names the file a maintainer would have to open, not the
        # header that happens to carry the USE
        catalog = {
            "a.F": cd.FileInfo(Path("a.F"), {"b.h"}, set(), set(), True),
            "b.h": cd.FileInfo(Path("b.h"), set(), {"gone"}, set(), False),
        }
        notes = cd.find_unresolved_uses(catalog, {})
        self.assertEqual(len(notes), 1)
        self.assertIn("gone", notes[0])
        self.assertIn("a.F", notes[0])

    def test_notes_do_not_change_the_emitted_dependencies(self):
        # find_unresolved_uses is pure diagnostics: an unresolvable
        # USE must not turn into a prerequisite, or make fails with
        # "No rule to make target" on a module it cannot build
        catalog = {
            "a.F": cd.FileInfo(Path("a.F"), set(), {"xios"}, set(), True),
        }
        self.assertEqual(cd.resolve_dependencies(catalog, {})["a.F"], set())


class TestWriteMakeDepend(unittest.TestCase):
    def test_output_shape(self):
        # the writer is deliberately dumb: one "target.o: source deps"
        # line per given file, plus the do-not-edit banner. By the
        # time deps reaches it, resolve_dependencies has already
        # converted any module dependency to .o form, so the writer
        # just emits what it is given, verbatim.
        out = Path(self.id() + ".depend")
        self.addCleanup(out.unlink)
        deps = {"a.F": {"b.h", "modc.o", "raw_included.F"}}
        cd.write_make_depend(out, [Path("a.F")], deps, [])
        text = out.read_text()
        self.assertIn("a.o: a.F", text)
        self.assertIn("b.h", text)
        self.assertIn("modc.o", text)
        self.assertIn("raw_included.F", text)
        self.assertIn("DO NOT EDIT", text)


class TestMainEndToEnd(unittest.TestCase):
    def test_no_args_errors(self):
        # `make depend` passes a *.F90 *.F glob; if it ever expands to
        # nothing, fail loudly rather than quietly writing an empty
        # Make.depend that disables every dependency in the build
        self.assertEqual(cd.main([]), 1)

    def test_help_flag_exits_zero_without_needing_files(self):
        # --help has to work before argument validation, so it is
        # usable from anywhere, without a source tree
        self.assertEqual(cd.main(["--help"]), 0)
        self.assertEqual(cd.main(["-h"]), 0)

    def test_fixture_tree(self):
        # the whole pipeline on a two-file tree: scan, register the
        # module, resolve, write. Also covers the default output path,
        # Make.depend in the current directory.
        tmp = Path(self.id())
        tmp.mkdir(exist_ok=True)
        self.addCleanup(lambda: [p.unlink() for p in tmp.glob("*")] and tmp.rmdir())
        (tmp / "modb.F90").write_text("      module modb\n      end module\n")
        (tmp / "a.F90").write_text("      use modb\n      Use modb\n")
        import os
        cwd = os.getcwd()
        os.chdir(tmp)
        try:
            self.assertEqual(cd.main(["a.F90", "modb.F90"]), 0)
            out = Path("Make.depend").read_text()
            self.assertIn("a.o:", out)
            self.assertIn("modb.o", out)
        finally:
            os.chdir(cwd)
            Path(tmp / "Make.depend").unlink(missing_ok=True)

    def test_output_option_writes_elsewhere(self):
        # -o keeps the tool out of the source tree: the Makefile names
        # the file explicitly, and check_differential.py can run both
        # tools without clobbering a real Make.depend
        tmp = Path(self.id())
        tmp.mkdir(exist_ok=True)
        self.addCleanup(lambda: [p.unlink() for p in tmp.glob("*")] and tmp.rmdir())
        (tmp / "modb.F90").write_text("      module modb\n      end module\n")
        (tmp / "a.F90").write_text("      use modb\n")
        import os
        cwd = os.getcwd()
        os.chdir(tmp)
        try:
            for flag in (["-o", "custom.d"], ["--output", "custom.d"],
                         ["--output=custom.d"]):
                Path("custom.d").unlink(missing_ok=True)
                self.assertEqual(cd.main([*flag, "a.F90", "modb.F90"]), 0)
                self.assertTrue(Path("custom.d").exists(), flag)
                self.assertIn("modb.o", Path("custom.d").read_text())
                self.assertFalse(Path("Make.depend").exists(), flag)
        finally:
            os.chdir(cwd)

    def test_output_option_without_argument_errors(self):
        # a dangling -o is a usage error, not an IndexError traceback
        self.assertEqual(cd.main(["-o"]), 1)

    def test_no_fast_filter_flag_produces_same_result(self):
        # the escape hatch is only worth keeping if the two paths are
        # interchangeable end to end, not just statement by statement
        # (TestFilterNeverStricterThanParser covers that half)
        tmp = Path(self.id())
        tmp.mkdir(exist_ok=True)
        self.addCleanup(lambda: [p.unlink() for p in tmp.glob("*")] and tmp.rmdir())
        (tmp / "modb.F90").write_text("      module modb\n      end module\n")
        (tmp / "a.F90").write_text("      use modb\n      Use modb\n")
        import os
        cwd = os.getcwd()
        os.chdir(tmp)
        try:
            self.assertEqual(cd.main(["a.F90", "modb.F90"]), 0)
            fast_out = Path("Make.depend").read_text()
            self.assertEqual(
                cd.main(["--no-fast-filter", "a.F90", "modb.F90"]), 0
            )
            full_out = Path("Make.depend").read_text()
            self.assertEqual(fast_out, full_out)
        finally:
            os.chdir(cwd)
            Path(tmp / "Make.depend").unlink(missing_ok=True)


if __name__ == "__main__":
    unittest.main()

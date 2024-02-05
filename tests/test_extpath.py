
from unittest import TestCase
from astroscripts import extpath
from pathlib import PurePosixPath
import os

class TestExtFileName(TestCase):
    def test_creation_onlyname(self):
        extname = extpath.ExtFileName('somefile.fts')
        self.assertEqual(str(extname), 'somefile.fts')

        extname = extpath.ExtFileName('somefile.fts', hdu=None, filter=None)
        self.assertEqual(str(extname), 'somefile.fts')

        with self.assertRaises(ValueError):
            extpath.ExtFileName(PurePosixPath('somefile.fts'), hdu=None, filter=None)

    def test_creation_hdu(self):
        extname = extpath.ExtFileName('somefile.fts', hdu=2)
        self.assertEqual(str(extname), "somefile.fts[2]")

        extname = extpath.ExtFileName('somefile.fts', hdu='GTI')
        self.assertEqual(str(extname), "somefile.fts['GTI']")

        with self.assertRaises(TypeError):
            extpath.ExtFileName('somefile.fts', hdu=2.45)

    def test_creation_filter(self):
        flt="""X=Y*4.34+'fgh'"""
        extname = extpath.ExtFileName('somefile.fts', filter=flt)
        self.assertEqual(str(extname), "somefile.fts[{!r}]".format(flt))

        flt="""X=Y/5"""
        extname = extpath.ExtFileName('somefile.fts', filter=flt)
        self.assertEqual(str(extname), "somefile.fts[{!r}]".format(flt))

    def test_repr(self):
        extname = extpath.ExtFileName('somefile.fts', hdu='GTI')
        self.assertEqual(repr(extname), "<ExtFileName('somefile.fts', hdu='GTI', filter=None)>")

    def test_from_string(self):
        extname = extpath.ExtFileName.from_string("""somefile.fts[GTI]['X=Y']""", without_hdu=False)
        self.assertEqual(str(extname), "somefile.fts['GTI']['X=Y']")

        extname = extpath.ExtFileName.from_string("""somefile.fts['GTI']['X="Y"*5']""", without_hdu=False)
        self.assertEqual(str(extname), """somefile.fts['GTI']['X="Y"*5']""")

        extname = extpath.ExtFileName.from_string("""somefile.fts['X="Y"*5']""", without_hdu=True)
        self.assertEqual(str(extname), """somefile.fts['X="Y"*5']""")

    # def test_slice(self):
    #     extname = extpath.ExtFileName('somefile.fts', hdu=2)
    #     self.assertEqual(str(extname[3:]), "efile.fts[2]")
    #     self.assertEqual(str(extname[:4]), "some[2]")


    # def test_str_add_self(self):
    #     extname = extpath.ExtFileName('somefile.fts', hdu=2)
    #     self.assertEqual(str('text_'+extname), "text_somefile.fts")
    #
    #     extname = extpath.ExtFileName('somefile.fts', hdu=2)
    #     self.assertEqual('text_'+str(extname), "text_somefile.fts[2]")
    #
    # def test_self_add_str(self):
    #     extname = extpath.ExtFileName('somefile.fts', hdu=2)
    #     self.assertEqual(str(extname+'_text'), "somefile.fts_text[2]")
    #
    #     extname = extpath.ExtFileName('somefile.fts', hdu=2)
    #     self.assertEqual(str(extname)+'_text', "somefile.fts[2]_text")


class TestExtPath(TestCase):
    def __init__(self, *args):
        super().__init__(*args)
        self.extname = extpath.ExtFileName('somefile.fts', hdu=2, filter='X="Y"*5')
    def _create(self):
        return extpath.ExtPath('/home/work', self.extname)

    def _create_rel(self):
        return extpath.ExtPath('./work', self.extname)
    def test_create(self):
        self.assertEqual(str(self._create_rel()), """work/somefile.fts[2]['X="Y"*5']""")
        self.assertEqual(str(self._create()), """/home/work/somefile.fts[2]['X="Y"*5']""")
        self.assertEqual(repr(self._create()), "ExtPathAbs('/home/work'/{!r})".format(self.extname))
        self.assertEqual(type(self._create_rel()), extpath.ExtPath)
        self.assertEqual(type(self._create()), extpath.ExtPathAbs)

    def test_create_fromstr(self):
        self.assertEqual(str(extpath.ExtPath("""/home/work/somefile.fts[2]['X="Y"*5']""")),
                         """/home/work/somefile.fts[2]['X="Y"*5']""")
        self.assertEqual(str(extpath.ExtPath("""/home/work/somefile.fts[2]['X="Y"*5']""", hdu=33)),
                         """/home/work/somefile.fts[33]""")

    def test_create_fromextpath(self):
        expa = self._create_rel()
        self.assertEqual(str(extpath.ExtPath('/bin', expa)), """/bin/work/somefile.fts[2]['X="Y"*5']""")
        self.assertEqual(str(extpath.ExtPath('./bin', expa)), """bin/work/somefile.fts[2]['X="Y"*5']""")
        self.assertEqual(type(extpath.ExtPath('./bin', expa)), extpath.ExtPath)
        self.assertEqual(type(extpath.ExtPath('/bin', expa)), extpath.ExtPathAbs)
    def test_create_fromfilepath(self):
        fpa = extpath.FilePath('work/somefile.fts')
        self.assertEqual(str(extpath.ExtPath('/bin', fpa, hdu=2)), """/bin/work/somefile.fts[2]""")
        self.assertEqual(str(extpath.ExtPath('./bin', fpa, hdu=2)), """bin/work/somefile.fts[2]""")
        self.assertEqual(str(extpath.ExtPath('./bin', fpa, filter='X=Y')), """bin/work/somefile.fts['X=Y']""")
        self.assertEqual(type(extpath.ExtPath('./bin',fpa)), extpath.ExtPath)
        self.assertEqual(type(extpath.ExtPath('/bin', fpa)), extpath.ExtPathAbs)

    def test_create_fromextname(self):
        self.assertEqual(str('./work' / self.extname), """work/somefile.fts[2]['X="Y"*5']""")
        self.assertEqual(type('./work' / self.extname), extpath.ExtPath)
        self.assertEqual(type('/work' / self.extname), extpath.ExtPathAbs)

    def test_fspath(self):
        self.assertEqual(os.fspath(self._create()), """/home/work/somefile.fts""")
    def test_name(self):
        self.assertEqual(str(self._create().name), """somefile.fts""")
    def test_extname(self):
        self.assertEqual(str(self._create().extname), """somefile.fts[2]['X="Y"*5']""")
    def test_stem(self):
        self.assertEqual(str(self._create().stem), """somefile""")
    def test_suffix(self):
        self.assertEqual(str(self._create().suffix), """.fts""")
    def test_parent(self):
        self.assertEqual(str(self._create().parent), """/home/work""")

    def test_withname(self):
        self.assertEqual(str(self._create_rel().with_name('newname.flc').name), """newname.flc""")
        self.assertEqual(str(self._create().with_name('newname.flc').name), """newname.flc""")
        self.assertEqual(type(self._create_rel().with_name('newname.flc')), extpath.ExtPath)
        self.assertEqual(type(self._create().with_name('newname.flc')), extpath.ExtPathAbs)
    def test_withextname(self):
        extname2 = extpath.ExtFileName('newname.gti', hdu='GTI')
        self.assertEqual(str(self._create().with_extname(extname2).extname),
                         """newname.gti['GTI']""")
        self.assertEqual(type(self._create_rel().with_extname(extname2)), extpath.ExtPath)
        self.assertEqual(type(self._create().with_extname(extname2)), extpath.ExtPathAbs)
    def test_withstem(self):
        self.assertEqual(str(self._create().with_stem('newname').name), """newname.fts""")
        self.assertEqual(type(self._create_rel().with_stem('newname')), extpath.ExtPath)
        self.assertEqual(type(self._create().with_stem('newname')), extpath.ExtPathAbs)
    def test_withsuffix(self):
        self.assertEqual(str(self._create().with_suffix('.ogg').name), """somefile.ogg""")
        self.assertEqual(type(self._create_rel().with_suffix('.ogg')), extpath.ExtPath)
        self.assertEqual(type(self._create().with_suffix('.ogg')), extpath.ExtPathAbs)
    def test_withparent(self):
        self.assertEqual(str(self._create().with_parent('/var/bin')), """/var/bin/somefile.fts[2]['X="Y"*5']""")
        self.assertEqual(type(self._create_rel().with_parent('/var/bin')), extpath.ExtPathAbs)
        self.assertEqual(type(self._create().with_parent('/var/bin')), extpath.ExtPathAbs)

    def test_truediv(self):
        self.assertEqual(str('/bin'/ self._create()), """/home/work/somefile.fts[2]['X="Y"*5']""")
        self.assertEqual(str('/bin'/ self._create_rel()), """/bin/work/somefile.fts[2]['X="Y"*5']""")
        self.assertEqual(type('./bin' / self._create()), extpath.ExtPathAbs)  # './bin' will be ignored
        self.assertEqual(type('./bin' / self._create_rel()), extpath.ExtPath)
        self.assertEqual(type('/bin' / self._create_rel()), extpath.ExtPathAbs)

cwd = os.getcwd() + '/'
class TestExtPathAbs(TestCase):
    def __init__(self, *args):
        super().__init__(*args)
        self.extname = extpath.ExtFileName('somefile.fts', hdu=2, filter='X="Y"*5')
    def _create(self):
        return extpath.ExtPath('./work', self.extname).absolute()
    def test_create(self):
        self.assertEqual(str(self._create()), cwd+"""work/somefile.fts[2]['X="Y"*5']""")
        self.assertEqual(repr(self._create()), "ExtPathAbs('{}work'/{!r})".format(cwd,self.extname))
    def test_cantcreate(self):
        with self.assertRaises(ValueError):
            extpath.ExtPathAbs('./work', self.extname)
    def test_fspath(self):
        self.assertEqual(os.fspath(self._create()), cwd+"""work/somefile.fts""")
    def test_name(self):
        self.assertEqual(str(self._create().name), """somefile.fts""")
    def test_extname(self):
        self.assertEqual(str(self._create().extname), """somefile.fts[2]['X="Y"*5']""")
    def test_stem(self):
        self.assertEqual(str(self._create().stem), """somefile""")
    def test_suffix(self):
        self.assertEqual(str(self._create().suffix), """.fts""")
    def test_parent(self):
        self.assertEqual(str(self._create().parent), cwd+"""work""")
        self.assertEqual(type(self._create().parent), extpath.FilePathAbs)

    def test_withname(self):
        self.assertEqual(str(self._create().with_name('newname.flc').name), """newname.flc""")
    def test_withextname(self):
        extname2 = extpath.ExtFileName('newname.gti', hdu='GTI')
        self.assertEqual(str(self._create().with_extname(extname2).extname),
                         """newname.gti['GTI']""")
    def test_withstem(self):
        self.assertEqual(str(self._create().with_stem('newname').name), """newname.fts""")
    def test_withsuffix(self):
        self.assertEqual(str(self._create().with_suffix('.ogg').name), """somefile.ogg""")
    def test_withparent(self):
        self.assertEqual(str(self._create().with_parent('/var/bin')), """/var/bin/somefile.fts[2]['X="Y"*5']""")
        with self.assertRaises(ValueError):
            self._create().with_parent('./var/bin')
    def test_withpap(self):
        self.assertEqual(str(self._create().with_parent_append('/var/bin')), '/var/bin'+cwd+"""work/somefile.fts[2]['X="Y"*5']""")
        with self.assertRaises(ValueError):
            self._create().with_parent('./var/bin')

    def test_truediv(self):
        #str must be ignored
        self.assertEqual(str('/bin' / self._create()), cwd+"""work/somefile.fts[2]['X="Y"*5']""")
        self.assertEqual(str('./bin' / self._create()), cwd+"""work/somefile.fts[2]['X="Y"*5']""")



    # def test_name(self):
    #         self.assertEqual(str(self._create().name), """somefile.fts""")

    # print(repr('/qww' / expa))
    # print(repr(expa / 'qww'))


# extname = extpath.ExtFileName('somefile.fts', hdu=2)
# extname2 = extpath.ExtFileName('qwefile.fts', hdu=2)
# print(str(extname))
# print(repr(extname))
# # print('%'.join(['text', extname]))
# expa = extpath.ExtPath('home/work', extname)
# print(expa)
# print(repr(expa))
# print(repr(expa.with_suffix('.ggt')))
# print(repr(expa.with_name(extname2)))

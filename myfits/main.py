"""General purpose functions."""

import os, re, copy
from typing import Union, Self, Callable, Type
from pathlib import PurePath
from astropy.io import fits
import mypython as my
from mypython import Actions, FilePath, FilePathAbs


class TaskError(Exception):
    """Exception to raise in custom scripts."""

    def __init__(self, taskname, custom_message=None, filename=''):
        self.taskname = taskname
        if custom_message:
            self.msg = custom_message
        else:
            self.msg = f"Task '{taskname}' caused an error"
            self.msg += f' during processing file {filename}.' if filename else '.'
        self.filename=filename
        super().__init__(self.msg)


#########################################################
#### Classes to provide paths with HDU and table columns


class ExtFileName:
    """Extended filename to store also hdu number and columns of FITS/ASCII tables.

    Notes. The PurePath's constructor calls str(os.fspath(obj)) for each
    of the input arguments and then uses the resultant string with '/'.join().
    The function os.fspath() do nothing for subclasses of str and calls
    __fspath__() otherwise. Therefore, if ExtFileName is a subclass of str,
    it must return the pure filename (i.e. without hdu or filter) via
    its __str__ method. If it is of any other type, the filepath must be
    returned by __fspath__, and the magic __str__ may have any output.
    """

    def __init__(self, name: str, hdu: Union[int, str] = None, filter=None):
        self._name = self._name_check(name)
        self.hdu = hdu
        self.filter = filter
        
    @classmethod
    def _name_check(cls, name):
        if not isinstance(name, str) or \
                my.FilePath(name).name != name:
            raise ValueError(f"Incorrect name '{repr(name)}'. Only pure "
                             "string names are allowed (not paths or other objects).")
        if '[' in name or ']' in name:
            my.printwarn(f"you are using square brackets in name: {name}.",
                         cls.__name__)
        return name

    def __str__(self):
        lst = [self._name]
        if self._hdu:
            lst.append("[{!r}]".format(self._hdu))
        if self.filter:
            lst.append("[{!r}]".format(self.filter))
        return ''.join(lst)

    def __repr__(self):
        return '<{}({!r}, hdu={!r}, filter={!r})>'.format(
            self.__class__.__name__, self._name, self._hdu, self.filter)

    def __fspath__(self):
        """Return filename for os.fspath()."""
        return self._name

    def __rtruediv__(self, other):
        if isinstance(other, (str, os.PathLike)):
            return ExtPath(other, self)
        return NotImplemented

    def get_opts(self):
        """Return hdu and filter arguments as a dict."""
        return dict(hdu=self._hdu, filter=self.filter)

    @property
    def name(self):
        return self._name

    @property
    def hdu(self):
        return self._hdu

    @hdu.setter
    def hdu(self, val):
        if not isinstance(val, (int, str, type(None))):
            raise TypeError("HUD key must be either a integer or a string")
        self._hdu = val

    # @property
    # def filter(self):
    #     return self._filter

    # @filter.setter
    # def filter(self, val):
    #     self._filter = val

    @classmethod
    def from_string(cls, expression: str, without_hdu=False):
        """Create the object from the parse_string() result."""
        return cls(*cls.parse_string(expression, without_hdu))

    @staticmethod
    def _process_filter(match1, without_hdu):
        if without_hdu is True:
            filter = match1.group(2)
            if match1.group(3):
                raise ValueError(f"Inadmissible expression {match1.group(3)} "
                                 "in the second pair of brackets.")
        else:
            filter = match1.group(3)

        if filter and len(filter) > 2 and filter[0] == filter[-1] and \
                (filter[0] == '"' or filter[0] == "'"):
            filter = filter[1:-1]

        return filter

    @classmethod
    def parse_string(cls, expression: str, without_hdu: bool = False,
                     only_path: bool = False):
        """Parse if expression has a format "file.fts[hdu]['filter']".

        Here hdu may be an integer representing an ordering number of the HDU in
        the FITS file or a string representing the HDU's name. The existence
        of the HDU won't be checked. Path and filter can be any strings.

        :param expression: String expression to analyse.
        :param without_hdu: If True it allows to omit the hdu key and
            consider the string in the first square brackets as
            the filtering exporession.
        :param only_path: Ignore the hdu and filter keys and return
            only the filepath. It may be useful for isolation of
            the pure filepath from the total expression without
            printing warnings.
        :return: The tuple of three elements: tuple(path, hdu, filter),
            where the path van be empty string while hdu and filter
            can be None.
        :raises ValueError:when the hdu string cannot be parsed.
        """
        hdu, filter = None, None
        match1 = re.match("(.*?)(?:\[(.*?)])?(?:\[(.*?)])?$", expression)  # It always gives some result,
        # at least match1.group(1) always exists but may be not a filename
        if only_path:
            return match1.group(1), hdu, filter

        if without_hdu is False:
            if match1.group(2):  # Parse first brackets
                hdustr = match1.group(2).replace('"', '').replace("'", '')  # remove all quotes
                if match2 := re.match("""^(\d+)|(\w+)$""", hdustr):  # int or string
                    if match2.group(1):
                        hdu = int(match2.group(1))
                    elif match2.group(2):
                        hdu = str(match2.group(2))
                else:  # Can't parse hdu key
                    without_hdu = True
                    my.printwarn(f"the string {match1.group(2)} will be "
                                 "interpreted as the filtering expression.", cls.__name__)

        filter = cls._process_filter(match1, without_hdu)
        return match1.group(1), hdu, filter


class ExtPath(FilePath):
    """Extend FilePath to store FITS's hdu and filtering expression

    The object can be constructed from path, path + ExtName, string or
    other ExtPath, where path is any os.Pathlike object; HDU and
     filter can be passed as keyword arguments."""

    def __new__(cls, *args: Union[str, PurePath, ExtFileName], **kwargs):
        superargs = list(args) if len(args)>0 else ['']
        last = superargs[-1]  # the last argument
        if isinstance(last, ExtFileName):
            extname = ExtFileName(last.name, **kwargs) if kwargs else copy.copy(last)
        elif isinstance(last, PurePath):
            extname = copy.copy(last._extname) if hasattr(last, '_extname') else \
                ExtFileName(last.name, **kwargs)
        elif isinstance(last, str):  # Try to parse string
            # Split to the pure string path and its extended part
            strpath, hdu, filter = ExtFileName.parse_string(last, only_path=bool(kwargs))
            extname_extra_args = kwargs or {'hdu':hdu, 'filter':filter}
            extname = ExtFileName(my.FilePath(strpath).name, **extname_extra_args)
            superargs[-1] = strpath   # save the string without the tail [hdu][filter]
        else:
            raise TypeError(f"Unsupported type of argument: {type(last)}.")
        obj = super().__new__(cls, *superargs)
        obj._extname = extname
        return obj

    def _from_parsed_parts(self, drv, root, parts):
        """Construct object from parts.

        This function overrides the PurePath's classmethod(!) used to
        create derivative subpaths from the original path (for self.with_name,
        self.parent, etc.). Now, it will be an ordinary method which gets extname
        from 'self' and puts it into the child object.
        """
        ancestor={ExtPath: FilePath, ExtPathAbs: FilePathAbs}
        cls = self.__class__
        name = os.fspath(parts[-1]) if len(parts) > 0 else ''
        #   for self.parent                   for self.with_parent
        if len(parts) != len(self._parts) and name != self.name:
            cls = ancestor[cls]
        obj = FilePath._from_parsed_parts.__func__(cls, drv, root, parts)
        obj._extname = self._extname
        return obj

    def __fspath__(self):
        """Return the normal path for os.fspath()."""
        return super().__str__()

    def __str__(self):
        """Return string representation."""
        tail = str(self._extname)[len(self._extname.name):]
        return super().__str__() + tail

    def __repr__(self):
        return "{}({!r}/{})".format(self.__class__.__name__,
                                    str(self.parent), repr(self._extname))

    def __truediv__(self, other):
        """Return self / other."""
        raise TypeError(f"Unsupported operation for {self.__class__}")

    def __rtruediv__(self, other):
        """Return other / self."""
        if isinstance(other, (str, os.PathLike)):
            return self.__class__(other, self)

    @property
    def extname(self):
        """Return name as a ExtFileName object."""
        return self._extname

    @property
    def hdu(self):
        return self._extname.hdu

    def with_extname(self: Self, new_extname: ExtFileName) -> Self:
        return self.__class__(self.parent, new_extname)


class ExtPathAbs(ExtPath, FilePathAbs):
    pass


ExtPath._class_abspath = ExtPathAbs
ExtPathAbs._class_abspath = ExtPathAbs


###############################################
####  Checks of the file format ###############


def _helper_check_file_is(func_deep_check: Callable, filetype: str,
                          filepath: Union[str, FilePath, ExtPath],
                          default_hdu: Union[str, int], action: Actions, progname) -> ExtPathAbs | None:
    extpath = ExtPath(filepath).absolute()
    if extpath.extname.hdu is None: extpath.extname.hdu = default_hdu
    hdu = extpath.extname.hdu
    try:
        with fits.open(extpath) as fts:
            if isinstance(hdu, int):    # hdu argument is a number (index)
                if hdu >= len(fts):
                    action.do(f"'{extpath.name}' doesn't have an extension #{hdu}.",
                                      exception_class=IndexError, progname=progname)
                    return None
            elif isinstance(hdu, str):  # hdu argument is a dict keyword
                if hdu not in fts:
                    action.do(f"'{extpath.name}' doesn't have an extension "
                              f"named '{hdu}'.", exception_class=KeyError, progname=progname)
                    return None
            else:
                raise TypeError("Wrong type of the 'hdu' argument")

            hduobj=fts[hdu]
            _deepcheck_action = Actions.WARNING if action is not Actions.NOTHING else Actions.NOTHING
            if func_deep_check(hduobj, action=_deepcheck_action):
                return extpath
            else:
                hdustr = f'#{hdu}' if isinstance(hdu, int) else f"'{hdu}'"
                action.do(f"'{extpath.name}' doesn't seem to be a valid {filetype} "
                          f"(extension {hdustr}).", exception_class=TypeError, progname=progname)
    except FileNotFoundError:
        action.do(f"File '{extpath.name}' is not found.", exception_class=FileNotFoundError,
                  progname=progname)
    except OSError:
        action.do(f"'{extpath.name}' is not a FITS file.", exception_class=TypeError,
                  progname=progname)
    return None


def fits_check_file_is_fits(filepath: Union[ExtPath, FilePath, str],
                            action: Union[Actions, str] = Actions.DIE, progname=None) -> ExtPathAbs | None:
    """Check if the file is readable FITS and return its absolute path."""
    return _helper_check_file_is(lambda x: True, 'FITS', filepath, 
                                 action=Actions(action), default_hdu=0, progname=progname)


def fits_check_file_is_lc(filepath: Union[ExtPath, FilePath, str],
                          action: Union[Actions, str] = Actions.DIE, progname=None) -> ExtPathAbs | None:
    """Check whether the first extension is a light curve."""
    return _helper_check_file_is(fits_check_hdu_is_lc, 'light curve', filepath,
                                 default_hdu=1, action=Actions(action), progname=progname)


def fits_check_file_is_gti(filepath: Union[ExtPath, FilePath, str],
                           action: Union[Actions, str] = Actions.DIE, progname=None) -> ExtPathAbs | None:
    """Check whether the first extension is a light curve."""
    return _helper_check_file_is(fits_check_hdu_is_gti, 'GTI file', filepath,
                                 default_hdu=1, action=Actions(action), progname=progname)


def _helper_check_hdu_is(HDU, class1_key:str, type_string:str, classmax:int, action:Actions, **kwargs):
    if 'HDUCLASS' not in HDU.header:
        action.do("There is no 'HDUCLASS' keyword in the HDU header. "
                 "Probably this FITS file does not conform OGIP standard.", **kwargs)
        return False, {}
    if HDU.header['HDUCLASS'].upper() != 'OGIP':
        action.do("'HDUCLASS' keyword is not 'OGIP'.", **kwargs)
        return False, {}
    classkeys={}
    for classkey in ['HDUCLAS'+str(i) for i in range(1,classmax+1)]:
        if classkey not in HDU.header:
            action.do(f"There is no '{classkey}' keyword in the HDU header. "
                     "Probably this FITS file does not conform OGIP standard.", **kwargs)
            return False, {}
        classkeys[classkey] = HDU.header[classkey].upper()
    if HDU.header['HDUCLAS1'].upper() != class1_key:
        action.do("The 'HDUCLAS1' keyword is not '{}'. "\
                 "Probably this is not a {} extension.".format(class1_key, type_string), **kwargs)
        return False, classkeys
    return True, classkeys


def fits_check_hdu_is_spectrum(HDU, withclasskeys: bool = False,
                               action: Union[Actions, str] = Actions.WARNING,
                               exception_class : Type[Exception] = TypeError, **kwargs):
    """Check if the HDU is a spectrum folowing the OGIP standard."""
    checkres, classkeys = _helper_check_hdu_is(HDU, 'SPECTRUM', 'spectral', 3, Actions(action),
                                               exception_class = exception_class, **kwargs)
    if not checkres:
        return (False, classkeys) if withclasskeys else False
    #TODO: other checks
    return (True, classkeys) if withclasskeys else True


def fits_check_hdu_is_lc(HDU, withclasskeys=False,
                         action: Union[Actions, str] = Actions.WARNING,
                         exception_class: Type[Exception] = TypeError, **kwargs):
    """Check if the HDU is a light curve folowing the OGIP standard."""
    checkres, classkeys = _helper_check_hdu_is(HDU, 'LIGHTCURVE', 'light curve', 3, Actions(action),
                                               exception_class = exception_class, **kwargs)
    if not checkres: 
        return (False, classkeys) if withclasskeys else False
    if 'TIME' not in HDU.columns.names:
        Actions(action).do("There is no column 'TIME'.", exception_class = exception_class, **kwargs)
        return (False, classkeys) if withclasskeys else False
    if 'RATE' not in HDU.columns.names:
        Actions(action).do("There is no column 'RATE'.", exception_class = exception_class, **kwargs)
        return (False, classkeys) if withclasskeys else False
    return (True, classkeys) if withclasskeys else True


def fits_check_hdu_is_gti(HDU, withclasskeys=False,
                         action: Union[Actions, str] = Actions.WARNING,
                         exception_class: Type[Exception] = TypeError, **kwargs):
    """Check if the HDU is a GTI-file folowing the OGIP standard."""
    checkres, classkeys = _helper_check_hdu_is(HDU, 'GTI', 'GTI', 2, Actions(action),
                                               exception_class = exception_class, **kwargs)
    if not checkres:
        return (False, classkeys) if withclasskeys else False
    if 'START' not in HDU.columns.names:
        Actions(action).do("There is no column 'START'.", exception_class = exception_class, **kwargs)
        return (False, classkeys) if withclasskeys else False
    if 'STOP' not in HDU.columns.names:
        Actions(action).do("There is no column 'STOP'.", exception_class = exception_class, **kwargs)
        return (False, classkeys) if withclasskeys else False
    if len(HDU.data) == 0:
        Actions(action).do("File seems to have a valid GTI extension, but it doesn't "
                          "contain any data.", exception_class = exception_class, **kwargs)
        return (False, classkeys) if withclasskeys else False
    return (True, classkeys) if withclasskeys else True


#### Fits tasks ######


def fits_keywords_getmany(HDU, keynames: list[str], defaults: dict = None, *,
                          action=Actions.EXCEPTION) -> dict[str]:
    """Read multiple keyword from the FITS header.

    The aim of this function to provide a functional analogues to dict.get()
    for headers' keywords.

    :param HDU: One object from the list obtained via astropy.io.fits.open().
    :param keynames: List of the keywords' names to read from FITS header
    :param defaults: If the keyword is missing use value from this dict.
        The default is None.
    :param action: Perform 'action' if the keyword is missing both in the
        header and in the defaults. The default is EXCEPTION
    :return: dict with keywords and their values
    """
    result={}
    for key in keynames:
        if key in HDU.header:  # Search for KEY or a KEYI+KEYF pair
            result[key] = HDU.header[key]
        elif (key+'I' in HDU.header) and (key+'F' in HDU.header):
            result[key] = HDU.header[key+'I'] + HDU.header[key+'F']
        else:   # Try to get it from user dict
            if key not in (defaults or {}):
                Actions(action).do(f"Required keyword '{key}' is not found",
                                   exception_class=KeyError)
            else:
                result[key] = defaults[key]
    return result


def gti_get_limits(gtifile: ExtPath):
    """Return boundaries of the GTI interval.

    :param gtifile: Path to the GTI FITS file
    :returns: a pair of float values
    """
    with fits.open(gtifile) as ftsgti:
        start=ftsgti[gtifile.hdu].data['START']
        stop=ftsgti[gtifile.hdu].data['STOP']
    return [start[0], stop[-1]]     # Get starttime from the first row and stoptime from the last row
    
#################################################   
#### DS9  ####################################

def ds9_send_commands(ds9, command_list):
  for line in command_list:
    ds9.set(line)


def ds9_open(filename='',title='',workdir='',colormap='bb', scale='log', zoom='',extra_commands=[]):
  import pyds9
  # import time
  
  ds9 = pyds9.DS9(wait=20)
  # time.sleep(3)
  
  commands=[]
  
  if workdir: commands.append('cd {}'.format(workdir))
  if filename: commands.append('file {}'.format(filename))
  if title: commands.append('title text {}'.format(title))
  if zoom: commands.append('zoom to {}'.format(''))
  commands.append('scale {}'.format(scale))
  commands.append('cmap {}'.format(colormap))
  
  ds9_send_commands(ds9,commands+extra_commands)
  return ds9
  
def ds9_load_region(ds9, filename, color='', group='',clear=False):

  if clear:
    #ds9.set('regions select all')
    ds9.set('region delete')
    
  if color or group:  #This demands to edit the region file content
    with open(filename,'r') as f:
      content=f.readlines()
    extra_text=' #'
    if color:
      extra_text+=' color={}'.format(color)
    if group:
      extra_text+=' group={}'.format(color)
      
    for i in range(3,len(content)):   #Ignore three first header lines 
      content[i]=content[i].replace('\n',extra_text+'\n')
    
    return ds9.set('regions',''.join(content))
    
  return ds9.set('regions load {}'.format(filename))
    
  
def ds9_load_image(ds9,filename,regfile='',extra_commands=[]):
  res=ds9.set('file {}'.format(filename))
  if regfile: ds9_load_region(ds9,regfile)
  if extra_commands: ds9_send_commands(ds9,extra_commands)
  return res

##########################################
#### XSELECT ###############################
def __xselect_command_generator(data, header=None):
  '''Generate extrAction commands. The 'data' is a list of dicts. Each dict contains
  keys 'prodtype', 'regfile', 'outfile' and also optionally 'chkimg',  
  'exta_filters', 'extra_arguments' (arguments for the xslect's "extract") and
  'extra_commands' (called after other commands) . The 'header' is the dict 
  containing 'session_name' and  'evtfiles'. If header is None, function 
  returns only core command (without 'read events' and 'exit' lines).
  The 'regfile' may be None to extract from the full FOV. '''
  
  cmd_start = cmd_core = cmd_end = ''   #Empty strings
  
  for entry in data:
    check_switching_arguments(entry['prodtype'],['spectrum','curve','events','image'])
    if entry['regfile']:
      cmd_core+='filter region {regfile}\n'.format(**entry)
    if 'extra_filters' in entry:
      for curfilter in entry['extra_filters']:
        cmd_core+='filter {filter_type} "{expression}"\n'.format(**curfilter)
    cmd_core+='extract {prodtype} {extra_arguments}\nsave {prodtype} {outfile}\n\n'.format(**entry) #additional '\n' for 'extract events'
    if 'chkimg' in entry:
      cmd_core+='extract image\nsave image {chkimg}\n'.format(**entry)
    if 'extra_commands' in entry:
      for line in entry['extra_commands']:
        cmd_core+=line+'\n'
  
    if header:
      ses_name = 'xsel{}'.format(os.getpid()) if 'session_name' not in header else header['session_name']
      evtfiles_txt=' '.join(header['evtfiles']) #Convert list of paths to one string with spaces
      cmd_start='{}\nread events\n./\n{}\nyes\n'.format(ses_name, evtfiles_txt)  
      cmd_end='exit\nn\n'

  return cmd_start+cmd_core+cmd_end
  
  
def __xselect_extract_product_with_background(product_type, evtfiles, objprodname, objregfile, bkgregfile, bkgprodname, outdir, clobber, make_checkimages, savelog, extra_extraction_arguments='', extra_filter_commands=[]):

  if type(evtfiles) is not list:   
    raise TypeError("'evtfiles' must be a python list")

  if outdir:
    outdir+='/'  
   
  #Check the presence of input files and get abs pathes
  if product_type == 'events' and objregfile is None:
    objreg_abspath=None
  else:
    objreg_abspath = check_file_exists(objregfile,action='exception')
  res={'object':check_file_not_exist_or_remove(outdir+objprodname, overwrite=clobber,action='exception') }    #Return absolute paths of the extracted files
  data_for_generator=[{'prodtype':product_type,'outfile':objprodname,'regfile':objreg_abspath,
  'extra_arguments':extra_extraction_arguments, 'extra_filters':extra_filter_commands}]   #Entries for object
  
  if bkgregfile:          #If the backgroud spectrum/lcurve is also needed
    bkgreg_abspath = check_file_exists(bkgregfile,action='exception')
    if not bkgprodname:   #Use defautl name
      bkgprodname = name_generator(objprodname,{'ends_with':'_bkg'})   #Add '_bkg' to the object filename
    res['background']=check_file_not_exist_or_remove(outdir+bkgprodname, overwrite=clobber,action='exception')
    data_for_generator[0]['extra_commands']=['clear region']    #Remove object's region before background extrAction
    data_for_generator+=[{'prodtype':product_type,'outfile':bkgprodname,'regfile':bkgreg_abspath,
    'extra_arguments':extra_extraction_arguments, 'extra_filters':extra_filter_commands}]    #Append entries for background
    
    
  #Go to outdir
  startWD=os.getcwd()
  if outdir:
    if not os.path.exists(outdir):
      os.mkdir(outdir)
    os.chdir(outdir)
  curWD=os.getcwd()

  if make_checkimages:    #Add check-image entries to data_for_generator if needed
    for entry in data_for_generator:
      entry['chkimg'] = name_generator(entry['outfile'],{'starts_with':'tmpimg_', 'add_extension':'fts'})
      check_file_not_exist_or_remove(entry['chkimg'], overwrite=True)    #Just remove if the file exists
      
    
  #Convert paths of evt-files to relative
  prepared_evtlist=[]
  for evtfile in evtfiles:
    evt_reldir = os.path.relpath(os.path.dirname(evtfile),curWD)+'/'
    prepared_evtlist.append(evt_reldir+os.path.basename(evtfile))
    
  header_for_generator= {'evtfiles':prepared_evtlist}  #Write new evtlist here to pass the command generator
  
  #All is ready. Perform the extraction!
  callftools('xselect',stdin=__xselect_command_generator(data_for_generator,header_for_generator))
  
  for entry in data_for_generator:  #Convert all check-images to png
    if 'chkimg' in entry:
      fitstoimg(entry['chkimg'])
      os.remove(entry['chkimg'])
  
  if savelog:   #Save log if needed
    os.rename(outdir+'xselect.log',outdir+'log_{}.log'.format(objprodname))

  os.chdir(startWD)
  return res
      
def xselect_extract_spectrum(evtfiles, objspecname, objregfile, bkgregfile='', bkgspecname='', outdir='', keybackfile='', keyrespfile='', keyancrfile='', group=None, make_checkimages=True, clobber=False, savelog=False):

  extracted_files = __xselect_extract_product_with_background('spectrum', evtfiles, objspecname, objregfile, bkgregfile, bkgspecname, outdir, clobber, make_checkimages, savelog)

  #Determine fits-header keword BACKFILE
  keybackfile_ = keybackfile
  if 'background' in extracted_files:   #Background spectrum has been extracted
    keybackfile_ = os.path.basename(extracted_files['background'])    #Replace keyword
    if keybackfile:
      printwarn("Parameter 'keybackfile' will be replaces with newly extracted background spectrum")
    
  #Write keywords
  if keybackfile_:
    fits.setval(extracted_files['object'], 'BACKFILE', value=keybackfile_, ext=1)   
  if keyrespfile:
    fits.setval(extracted_files['object'], 'RESPFILE', value=keyrespfile, ext=1)   
  if keyancrfile:
    fits.setval(extracted_files['object'], 'ANCRFILE', value=keyancrfile, ext=1)   
    
  if group:
    extracted_files['grouped'] = group_spectrum_with_grppha(extracted_files['object'], group_min_counts=group, clobber=clobber)
    
  return extracted_files
  
def xselect_extract_lcurve(evtfiles, objlcurvename, objregfile, bkgregfile='', bkglcurvename='', outdir='', binsize=None, fracexp_thresh=0.95, pi_filter=[], make_checkimages=True, clobber=False, savelog=False):
  
  extra_arguments='exposure={} '.format(fracexp_thresh)
  if binsize:
    extra_arguments+='bin={}'.format(binsize)
    
  filter_commands=[]
  if pi_filter:
    filter_commands.append({'filter_type':'column', 'expression':'PI={:d}:{:d}'.format(pi_filter[0],pi_filter[1])})
    
  extracted_files = __xselect_extract_product_with_background('curve', evtfiles, objlcurvename, objregfile, bkgregfile, bkglcurvename, outdir, clobber, make_checkimages, savelog, extra_extraction_arguments=extra_arguments, extra_filter_commands=filter_commands)
  return extracted_files
  
def xselect_combine_evt(inevtfiles, outevtfile, outdir='', make_checkimages=True, clobber=False, savelog=False):
  extracted_files = __xselect_extract_product_with_background('events', inevtfiles, outevtfile, None, None, None, outdir, clobber, make_checkimages, savelog)
  return extracted_files
  
def xselect_filter_evt(inevtfiles, outevtfile, regfile, outdir='', make_checkimages=True, clobber=False, savelog=False):
  extracted_files = __xselect_extract_product_with_background('events', inevtfiles, outevtfile, regfile, None, None, outdir, clobber, make_checkimages, savelog)
  return extracted_files

  
  
########################################################################
def group_spectrum_with_grppha(infile, outfile='', group_min_counts=1, clobber=False):
  
  if not outfile:
    outfile=name_generator(infile,{'ends_with':'_grp{}'.format(group_min_counts)})   #Add _grpN to the infile name
  
  check_file_not_exist_or_remove(outfile, overwrite=clobber,action='exception')  
  if callftools('grppha {} {}'.format(infile,outfile), stdin='group min {}\nexit\n'.format(group_min_counts)):
    return outfile

  return None


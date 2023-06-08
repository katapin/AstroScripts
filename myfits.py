"""Manipulation with FITS files, general purpose functions."""

import sys, os, re
from contextlib import contextmanager
from astropy.io import fits
import mypython as my
from mypython import *

        

###################################################    
#### Callers of the external programms  ####################

        
def callshell(cmd, stdin='', separate_logfile='', return_code=False):
    """Call shell programm, pass input and log result."""
    return my.callandlog(cmd, stdin=stdin, separate_logfile=separate_logfile, 
                         return_code=return_code, progname='callshell')
  
def callftools(cmd, stdin='', separate_logfile='', return_code=False):
    """Call ftools programm, pass input and log result."""
    return my.callandlog(cmd, stdin=stdin, separate_logfile=separate_logfile,
           extra_start='source ~/.bashrc \nheainit\n', return_code=return_code, progname='callftools')

def callciao(cmd, stdin='', separate_logfile='', return_code=False):
    """Call Chandra ciao programm, pass input and log result."""
    return my.callandlog(cmd, stdin=stdin, separate_logfile=separate_logfile, 
          extra_start='source ~/.bashrc \nciaoinit >/dev/null\n', 
          return_code=return_code,  progname='callciao')

def logging_turn_on(logfile:str, progname:str=''):
    """Turn on the logginng in the user's scripts.
    
    Helper function that checks presence of the logfile, enables logging and
    writes starting message with sys.argv parameters of the script. If the file
    at provided path exists, it prints error and calls sys.exit().

    Parameters
    ----------
    logfile : str
        Path to logfile
    progname : str, optional
        Name of script to write it in the logfile. The default is ''.

    Returns
    -------
    None.

    """
    my.check_file_not_exist_or_remove(logfile,overwrite=False,action=my.Action.DIE,
             extra_text='Please use another name or remove it manually.')
    my.set_logfile(logfile)
    args = sys.argv.copy()
    args[0] = os.path.basename(args[0])
    my.logtext('Started as "{}"'.format(' '.join(args)), progname)
  
    
###############################################
####  Checks of the file format ##############
        
    
def _helper_check_file_is(func, filetype, filepath, hdu, **kwargs):
    try:
        with fits.open(filepath) as fts:
            if isinstance(hdu, int):
                if hdu >= len(fts):
                    my._do_action(text=f"'{filepath}' doesn't have an extension #{hdu}.", 
                          exception_class=IndexError, **kwargs)
            elif isinstance(hdu, str):
                if not hdu in fts:
                    my._do_action(text=f"'{filepath}' doesn't have an extension named '{hdu}'.", **kwargs)
            else:
                raise TypeError("Wrong type of the 'hdu' argument")

            hduobj=fts[hdu]
            if func(hduobj):
                return my.FilePath(filepath).make_absolute()
            else:
                my._do_action(text=f"'{filepath}' doesn't seem to be a valid {filetype} (extension #{fts.index(hduobj)}).", **kwargs)
    except FileNotFoundError:
        my._do_action(text=f"File '{filepath}' is not found.", exception_class=FileNotFoundError, **kwargs)
    except OSError:
        my._do_action(text=f"'{filepath}' is not a FITS file.", exception_class=OSError, **kwargs)

    return None
    
def fits_check_file_is_fits(filepath, action=my.Action.DIE, progname=None):
    """Check if the file is readable FITS and return its absolute path."""
    return _helper_check_file_is(lambda x: True, 'FITS', filepath, 
                                 action=action, hdu=0, progname=progname)

def fits_check_file_is_lc(filepath, hdu=1, action=my.Action.DIE, progname=None):
    """Check whether the first extension is a light curve."""
    return _helper_check_file_is(fits_check_hdu_is_lc, 'light curve', filepath, hdu=hdu, action=action, progname=progname)

def fits_check_file_is_gti(filepath, hdu=1, action=my.Action.DIE, progname=None):
    """Check whether the first extension is a light curve."""
    return _helper_check_file_is(fits_check_hdu_is_gti, 'GTI file', filepath, hdu=hdu, action=action, progname=progname)

def fits_check_hdu_is_lc(HDU):
    """Check if the HDU is a light curve folowing the OGIP standard."""
    if 'HDUCLASS' not in HDU.header:
        my.printwarn("There is no 'HDUCLASS' keyword in the HDU header. "
        "Probably this FITS file does not conform OGIP standard.")
        return False
    if HDU.header['HDUCLASS']!='OGIP':
        my.printwarn("'HDUCLASS' keyword is not 'OGIP'. "
        "Only OGIP FITS files supported.")
        return False
    if 'HDUCLAS1' not in HDU.header:
        my.printwarn("There is no 'HDUCLAS1' keyword in the HDU header. "
        "Probably this FITS file does not conform OGIP standard. ")
        return False
    if HDU.header['HDUCLAS1']!='LIGHTCURVE':
        my.printwarn("The 'HDUCLAS1' keyword is not 'LIGHTCURVE'. "
        "Probably this is not a light curve extension.")
        return False
    if not 'TIME' in HDU.columns.names:
        my.printwarn("There is no column 'TIME'.")
        return False
    if not 'RATE' in HDU.columns.names:
        my.printwarn("There is no column 'RATE'.")
        return False
    return True
    
def fits_check_hdu_is_gti(HDU):
    """Check if the HDU is a GTI-file folowing the OGIP standard."""
    if 'HDUCLASS' not in HDU.header:
        my.printwarn("There is no 'HDUCLASS' keyword in the HDU header. "
        "Probably this FITS file does not conform OGIP standard.")
        return False
    if HDU.header['HDUCLASS']!='OGIP':
        my.printwarn("'HDUCLASS' keyword is not 'OGIP'.")
        return False
    if 'HDUCLAS1' not in HDU.header:
        my.printwarn("There is no 'HDUCLAS1' keyword in the HDU header. "
        "Probably this FITS file does not conform OGIP standard.")
        return False
    if HDU.header['HDUCLAS1']!='GTI':
        my.printwarn("The 'HDUCLAS1' keyword is not 'GTI'. "
        "Probably this is not a GTI extension.")
        return False
    if not 'START' in HDU.columns.names:
        my.printwarn("There is no column 'START'.")
        return False
    if not 'STOP' in HDU.columns.names:
        my.printwarn("There is no column 'STOP'.")
        return False
    if len(HDU.data) == 0:
        my.printwarn("File seems to have a valid GTI extension, but it doesn't "
                     "contain any data")
        return False
    return True

def fits_parse_expression(file_and_hdu:str):
    """Split the FITS HDU from a filename given in the format 'file.fts[hdu]'.
    
    Here hdu may be an integer representing a ordering number of the HDU or 
    a string representing the HDU's name. If the expression can't be matched,
    it returns None. The existence of the hdu won't be checked.
    
    Parameters
    ----------
    file_and_hdu : str
        Expresion to analyse.

    Returns
    -------
    filepath : str
        DESCRIPTION.
    hdu : TYPE
        DESCRIPTION.

    """
    hdu=None
    match1=re.match("(.*?)(\[.*?\])?$", file_and_hdu) #It always gives some result,
    #at least match1.group(1) always exists may might be not a filename
    if match1.group(2): #This is the part with square brackets
        if match2:=re.match("^\[(\d+)\]?$", match1.group(2)):
            hdu = int(match2.group(1)) 
        elif match2:=re.match("^\['(\w+)'\]?$", match1.group(2)):            
            hdu = str(match2.group(1))
        elif match2:=re.match('^\["(\w+)"\]?$', match1.group(2)):            
            hdu = str(match2.group(1))
        else:
            return None
    filepath=match1.group(1)
    return filepath, hdu

def gti_get_limits(gtifile:str, hdu=1):
    """Return boundaries of the GTI interval.

    Parameters
    ----------
    gtifile : str
        Path to the GTI FITS file
    hdu : TYPE, optional
        Number of the HDU extension that contains GTI data. The default is 1.

    Returns
    -------
    typle
        Returns a pair of float values
    """
    with fits.open(gtifile) as ftsgti:
        start=ftsgti[hdu].data['START']
        stop=ftsgti[hdu].data['STOP']
    return [start[0], stop[-1]] #Star from the first row and stop from the last row
    
def fitslc_plot(lcraw:str, gtifile_with_hdu:str, outepsimg:str, lcbkg:str=None, bkgratio:float=None, 
                lcnet:str=None, binsize:int=500):
    """Make a eps image of the light curve using lcurve task of xronos.

    Parameters
    ----------
    lcraw : str
        Path to the raw object light curve (extracted from the source aperture
        but not background subtracted).
    gtifile_with_hdu : str
        GTI file to obtain X-axis limits for the plot. Maight in a form of gtifile.fts[#hdu]
    outepsimg : str
        Path to eps file to store the result
    lcbkg : str, optional
        Path to the light curve extracted from the background aperture. The default is None.
    bkgratio : float, optional
        Area correction factor for the background light curve. It's mandatory if
        the 'lcbkg' argument is provided. The default is None.
    lcnet : str, optional
        Path to the net (background subtracted) light curve. The default is None.
    binsize : int, optional
        Desired temporal resolution in seconds. The default is 500 sec.

    Returns
    -------
    bool
        Return True.

    """
    _ownname=my.getownname()
    outepsimg = my.FilePath(outepsimg).make_absolute()
    add_quotes = lambda s: "'" + s + "'"
    gtifile, hdu = fits_parse_expression(gtifile_with_hdu)
    if not hdu: hdu=1
    gtilimits=gti_get_limits(gtifile, hdu)
    nbin=(gtilimits[1]-gtilimits[0])/binsize+10 #Bins per interval 
    
    
    lcurves = [add_quotes(lcraw)]  #list which will help to generate the argument\
        #string for the lcurve task below
    
    #If lcbkg is needed, get bkgratio and calc new background light curve
    if lcbkg:
        if not bkgratio:
            raise TypeError("Argument 'lcbkg' requires 'bkgratio' but NoneType is given.")
        my.printgreen(f'BKGRATIO={bkgratio}')    
        tmplcbkg=my.TempFile.generate_from('lcbkg.fts')
        if not (callftools("ftcalc '{}' '{}' RATE RATE*{:f}".format(lcbkg, tmplcbkg, bkgratio)) or 
                callftools("ftcalc '{}' '{}' ERROR ERROR*{:f} clobber=yes ".format(tmplcbkg, tmplcbkg,
                bkgratio))):
            my.printerr("Cannot apply correction for the BACKSCALE")
            raise my.ExternalTaskError('ftcalc', filename=lcbkg, caller=_ownname)
        lcurves.append(add_quotes(tmplcbkg))
        
        
    #Convert GTI to XRONOS format
    tmpgtixron=my.TempFile.generate_from('tmpgti.wi')
    if not callftools(f"gti2xronwin -i '{gtifile}' -o '{tmpgtixron}'"):
        my.printerr("Can't convert GTI to XRONOS format")
        raise my.ExternalTaskError('gti2xronwin', filename=lcbkg, caller=_ownname)
    
    #Prepare PCO file
    tmppco=my.TempFile.generate_from('tmppco.pco')
    with open(tmppco, 'w') as fpco:
        fpco.write("plot off\n")
        fpco.write("Col 1 on 2\n")
        if lcbkg: fpco.write("Col 2 on 3\n")
        if lcnet: 
            fpco.write("Col 4 on 4\n")
            lcurves.append(add_quotes(lcnet))
        fpco.write("Lab T {}\n".format(outepsimg.basename))
        fpco.write("Lab Y Rate, counts/sec\n")
        fpco.write("R\n")
        fpco.write(f"H {outepsimg}/CPS\n")
    fpco.close()
    
    cmd="lcurve {:d} {} window='{}' dtnb={} nbint={} outfile=' ' plot=yes "\
    "plotdev='/dev/null/PS' plotdnum=3 plotfile='{}' </dev/null "\
    ">/dev/null".format(len(lcurves), ' '.join(lcurves), tmpgtixron, binsize, 
                        int(nbin), tmppco)
    if not callftools(cmd):
        my.printerr("Can't plot the light curve")
        raise my.ExternalTaskError('lcurve', filename=lcraw, caller=_ownname)
        
    if not outepsimg.file_exists():
        my.printerr(f"Something is going wrong: '{outepsimg.basename}' has not been created.", _ownname)
        raise my.TaskError(_ownname, filename=outepsimg)

    return True

    
def fitsimg_to_png(ftspath:str, imgpath:str='', ds9_extra_commands:str='', 
                   convert_extra_commands:str=''):
    """Convert FITS image to png with ds9."""
    _ownname=my.getownname()
    tmpimg  = my.TempFile.generate('.ps')
    imgpath = my.FilePath(imgpath) if imgpath else my.FilePath(ftspath).replace_extension('png')

    if not callshell("ds9 {} -nopanner -nomagnifier -noinfo -colorbar no -view buttons no -geometry 1360x768 "\
    "-zoom to fit -cmap bb -scale log exp 10000 -scale log {} -print destination file "\
    "-print filename {} -print -exit".format(ftspath,ds9_extra_commands,tmpimg)):
        my.printerr("Cannot call 'ds9'.",_ownname),
        raise my.ExternalTaskError('ds9', filename=ftspath, caller=_ownname)
           
    if not callshell("convert -density 300 {} {} -background white "\
    "-flatten -fuzz 5%% -trim {}".format(tmpimg, convert_extra_commands, imgpath)):
        my.printerr(f"Cannot call 'convert', '{imgpath.basename}' is not created.",_ownname)
        raise my.TaskError('convert', filename=imgpath.basename, caller=_ownname)
    
    if not os.path.exists(imgpath):
        my.printerr(f"Something is going wrong: '{imgpath.basename}' has not been created.", _ownname)
        raise my.TaskError(_ownname, filename=imgpath)
        
    return True
                
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

######################################  
#### XSPEC ###########################

class Xspec():
    __instance = None
    def __new__(cls, *args, **kwargs):
        """Make the class singleton."""
        if cls.__instance:
            return cls.__instance
        else:
            cls.__instance = super().__new__(cls)
            return cls.__instance
    
    def __init__(self, set_rebin=(5,5), statMethod='cstat', show_background=False, xaxis='keV', 
                 query='yes', nTries=5):
      import xspec
      self.__dict__.update(xspec.__dict__)
      with self.__silence():
          self.Plot.xAxis = xaxis
          self.Plot.background=show_background
          self.Plot.setRebin(*set_rebin)
          self.Fit.query=query
          self.Fit.statMethod = statMethod
          self.nTries=nTries
  
    def load_spectra(self, spectra: list , clear=True, quiet=False):
        """Load multiple spectra."""
        #TODO: relative paths of backfile, respfile etc. may cause conflicts with change_dir()
        startWD=os.getcwd()
        startChatter=self.Xset.chatter
        
        if clear:
          self.AllData.clear()
          
        res=[]  #Return list of xspec.Spectrum() objects
        for entry in spectra:
            os.chdir(os.path.dirname(entry.specfile))  #Go the spectrum host dir to load correct bkg and resp files from its FITS-header
            with self.__silence():  #Hide the output because it shows message with no bkg and resp files which may be confusing
                sp = self.Spectrum(entry.specfile)
                if entry.bkgfile:
                    sp.background=entry.bkgfile
                if entry.respfile:
                    sp.response=entry.respfile
                if entry.arffile:
                    sp.response.arf=entry.arffile
                if entry.ignore:
                    for expression in entry.ignore:
                        sp.ignore(expression)
                
                os.chdir(startWD)
            res.append(sp)
            if not quiet:
                sp.show()
    
        return res
    
    def standard_analysis(self, quiet=False):
        """Fit the model and return parameter values with errors for single spectrum."""
        free_parameters=[]
        for i in range(1,6):
            try:
                mo=self.AllModels(i)
                for k in range(1,mo.nParameters+1):
                    par=mo(k)
                    if not par.frozen: free_parameters.append(k)
            except:
                break
        
        curDelta=100
        curTry=0
        with self.__silence(quiet):
            while (curDelta >= self.Fit.criticalDelta) and (curTry < self.nTries):
                self.Fit.perform()
                first_fit = self.Fit.statistic
                self.Fit.error('1.0 {}'.format(' '.join(map(str,free_parameters))))
                self.Fit.perform()
                second_fit = self.Fit.statistic
                curDelta = first_fit - second_fit
                curTry+=1
            else:
                if curTry == self.nTries:
                    print("Can't fit the model. Tries are exhausted")
                    return None
  
    def load_single_spectrum(self, specfile:str, bkgfile:str='', respfile:str='', arffile:str='', ignore:list=[], clear=True, quiet=False):
      """Append new spectrum to already loaded."""
      spec = self.XspecSpecDescription(specfile, bkgfile, respfile, arffile, ignore)
      return self.load_spectra([spec], clear, quiet)[0]
  
    def backupModel(self, index:int=1):
        """Get ModelDescription object."""
        return self.XspecModelDescription(self.AllModels(index))
    
    def restoreModel(self, model_description):
        """Load model from the ModelDescription object."""
        startChatter=self.Xset.chatter
        self.Xset.chatter=0
        mdl=self.Model(model_description.expression)
        for i in range(model_description.nParameters):
            mdl(i+1).values = model_description[i].values
        
        self.Xset.chatter=startChatter
    
    @contextmanager
    def __silence(self, quiet=True):
        if quiet:
            startChatter = self.Xset.chatter
            self.Xset.chatter=0
            yield
            self.Xset.chatter=startChatter
        else:
            yield

    class XspecSpecDescription():
        def __init__(self, specfile:str, bkgfile:str='', respfile:str='', arffile:str='', ignore:list=[]):
            self.specfile = specfile
            self.bkgfile  = bkgfile
            self.respfile = respfile
            self.arffile  = arffile
            self.ignore   = ignore
        
    class XspecModelDescription():
        def __init__(self, model_object):
            if model_object:
                self.expression=model_object.expression
                self.nParameters=model_object.nParameters
                self.__parameters=[]
                self.__components={}
                for ind,comp_name in enumerate(model_object.componentNames):
                    comp = self.__Component(model_object.__getattribute__(comp_name), ind)
                    self.__components[comp_name] = comp
                    self.__parameters += comp.parameters
                    
        def __getitem__(self, key):
            if type(key) is int:
                return self.__parameters[key]
            elif type(key) is str:
                return self.__components[key]
            else:
                raise TypeError('key must be string or integer')
        
        def __repr__(self):
            return "<ModelDescription '{}' at {}>".format(self.expression, hex(id(self)))
        
        # def __del_comp(self, index):
        #     for par in self.__parameters:
        #         if par.comp_index == index:
        #             self.__parameters.remove(par)
            
        # def del_comp(self, index):
        #     new_instance = copy.deepcopy(self)
        #     new_instance.__del_comp(index)
        #     return new_instance
        
        def print_summary(self):
            print(f"Model '{self.expression}'")
            for par in self.__parameters:
                print(' {} {}:\t{}\t=\t{}'.format(par.comp.index, par.comp.name, par.name, str(par.values)))
                
        class __Component():
            def __init__(self, comp_object, index: int):
                self.__name=comp_object.name
                self.__index=index
                self.__parameters={}
                for par_name in comp_object.parameterNames:
                    self.__parameters[par_name]=self.__Parameter(comp_object.__getattribute__(par_name), self)
                    
            def __getitem__(self, key):
                if not type(key) is str:
                    raise TypeError('index must be string')
                return self.__parameters[key]
            
            def __repr__(self):
                return "<Component '{}' at {}>".format(self.name, hex(id(self)))
            
            @property
            def name(self):
                return self.__name
            
            @property
            def index(self):
                return self.__index
                
            @property
            def parameters(self):
                return list(self.__parameters.values())
    
            @property
            def parameter_names(self):
                return list(self.__parameters.keys())

            
            class __Parameter():
                def __init__(self,param_object, comp):
                    self.__name=param_object.name
                    self.__index=param_object.index
                    self.__comp = comp
                    self.values=param_object.values
                    
                def __repr__(self):
                    return "<Parameter {}={} of '{}' at {}>".format(self.name, self.values[0],
                        self.__comp.name, hex(id(self)))
                
                @property
                def name(self):
                    return self.__name
                
                @property 
                def index(self):
                    return self.__index
                
                @property 
                def comp(self):
                    return self.__comp 
        
  
  
########################################################################
def group_spectrum_with_grppha(infile, outfile='', group_min_counts=1, clobber=False):
  
  if not outfile:
    outfile=name_generator(infile,{'ends_with':'_grp{}'.format(group_min_counts)})   #Add _grpN to the infile name
  
  check_file_not_exist_or_remove(outfile, overwrite=clobber,action='exception')  
  if callftools('grppha {} {}'.format(infile,outfile), stdin='group min {}\nexit\n'.format(group_min_counts)):
    return outfile

  return None


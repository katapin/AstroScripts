#A few string about this program


import os
import subprocess
from astropy.io import fits
from mypython import * 


def fitsopen(filename):
    if not os.path.isfile(filename):
        die("File '%s' is not found" % filename)
    try:
        f = fits.open(filename)
    except OSError:
        die("'%s' is not a FITS file" % filename)
    return f
    
def fitscheck(filename):
    f = fitsopen(filename)
    f.close()
    
    
def callandlog(cmd,logfile='',stdin=''):
  printbold(cmd,'callftools')
  if logfile:
    child=subprocess.Popen(cmd+' | tee -a '+logfile,shell=True,stdin=subprocess.PIPE)
  else:
    child=subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE)
    
  while child.poll() is None:
   child.stdin.write(stdin.encode('ascii'))   #Convert string to binary
   child.stdin.flush()                        #Push the text
   child.wait()                               #and wait the result
   
  return not child.returncode
        
  
def callftools(cmd,logfile='',stdin=''):
  printbold(cmd,'callftools')
  
  if logfile:
    child=subprocess.Popen('source ~/.bashrc \nheainit\n' +
    cmd+' | tee -a '+logfile,shell=True,stdin=subprocess.PIPE)
  else:
    child=subprocess.Popen('source ~/.bashrc \nheainit\n' +
    cmd,shell=True,stdin=subprocess.PIPE)
    
  while child.poll() is None:
   child.stdin.write(stdin.encode('ascii'))   #Convert string to binary
   child.stdin.flush()                        #Push the text
   child.wait()                               #and wait the result
   
  return not child.returncode
    
    
def callciao(cmd,logfile=""):
    
    printbold(cmd,"callciao")
    if logfile:
        res=subprocess.call("source ~/.bashrc \nciaoinit >/dev/null\n" +
        cmd+" | tee -a "+logfile,shell=True)
    else:
        res=subprocess.call("source ~/.bashrc \nciaoinit >/dev/null\n" +
        cmd,shell=True)
    return not res

def getgtilimits(gtifile,hdu=1):
    ftsgti=fitsopen(gtifile)
    start=ftsgti[hdu].data['START']
    stop=ftsgti[hdu].data['STOP']
    ftsgti.close()
    return [start[0], stop[-1]]
    
def lcplot(lcobj,imglc,lcbkg,lcnet,gtifile,bkgratio,dtnb=500):
    
    gtilimits=getgtilimits(gtifile)
    if gtilimits:
        nbin=(gtilimits[1]-gtilimits[0])/dtnb+10 #Bins per interval 
    else:
        printerr("Corrupted GTI file","lcplot")
        return False
    
    #Get bkgratio and calc new background light curve
    tmplcbkg=".lcbkg.fts"
    if not (callftools("ftcalc '%s' '%s' RATE RATE*%s" % (lcbkg, 
    tmplcbkg, bkgratio))):
        printerr("Can not correct background for BACKSCALE","lcplot")
        return False
    callftools("ftcalc '%s' '%s' ERROR ERROR*%s clobber=yes " % (tmplcbkg,
    tmplcbkg, bkgratio))
        
    #Convert GTI to XRONOS format
    tmpgtixron=".tmpgti.wi"
    cmd="gti2xronwin -i '%s' -o '%s'" % (gtifile, tmpgtixron)
    if not callftools(cmd):
        printerr("Can't convert GTI to XRONOS format","lcplot")
        return False
    
    #Prepare PCO file
    tmppco=".tmppco.pco"
    fpco= open(tmppco, 'w')
    fpco.write("plot off\n")
    fpco.write("Col 1 on 2\n")
    fpco.write("Col 2 on 3\n")
    fpco.write("col 4 on 4\n")
    fpco.write("Lab T %s\n" % imglc)
    fpco.write("Lab Y Rate, counts/sec\n")
    fpco.write("R\n")
    fpco.write("H %s/CPS\n" % imglc)
    fpco.close()
    
    cmd="lcurve 3 %s %s %s window='%s' dtnb=%s nbint=%s outfile=' ' " \
    "plot=yes plotdev='/dev/null/PS' plotdnum=3 plotfile='%s' </dev/null " \
    ">/dev/null" % (lcobj, tmplcbkg, lcnet, tmpgtixron, dtnb, int(nbin), 
    tmppco)
    res=callftools(cmd)
    os.remove(tmplcbkg)
    os.remove(tmpgtixron)
    os.remove(tmppco)
    
    printgreen("BKGRATIO=%f" % bkgratio,"lcplot")    
    if not res:
        printerr("Can not plot light curves","lcplot")
        return False
    if not os.path.isfile(imglc):
        printerr("Something is going wrong: %s is not "
        "created." % imglc,"lcplot")
        return False
        
    printbold("Saved "+imglc,"lcplot")
    return True

    
def fitstoimg(fts,img='',ds9_extra_commands="",conver_extra_commands=""):
    tmpimg=".%d.ps" % os.getpid()
    if not img:
      img=name_generator(fts,{'replace_extension':'png'})
      
    if os.system("ds9 %s -nopanner -nomagnifier -noinfo -colorbar no -view buttons no -geometry 1360x768 "
    "-zoom to fit -cmap bb -scale log exp 10000 -scale log %s -print destination file " \
        "-print filename %s -print -exit" % (fts,ds9_extra_commands,tmpimg)):
        printerr("Cannot call 'ds9'.","fitstoimg")
        return False
           
    res=os.system("convert -density 300 %s %s -background white " \
    "-flatten -fuzz 5%% -trim %s" % (tmpimg, conver_extra_commands, img))
    os.remove(tmpimg)
    if res:
        printerr("Cannot call 'convert' (imagemagic), '%s' is " 
        "not created." % img,"fitstoimg")
        return False
    return True
    
def chkkey(HDU,key):
    return (key in HDU.header)

        
def chkislc(HDU):
    if not chkkey(HDU,'HDUCLASS'):
        printwarn("There is no 'HDUCLASS' keyword in the HDU header. "
        "Probably this FITS file does not conform OGIP standard. "
        "Only OGIP FITS files supported.","chkislc")
        return False
    if HDU.header['HDUCLASS']!='OGIP':
        printwarn("'HDUCLASS' keyword is not 'OGIP'. "
        "Only OGIP FITS files supported.","chkislc")
        return False
    if not chkkey(HDU,'HDUCLAS1'):
        printwarn("There is no 'HDUCLAS1' keyword in the HDU header. "
        "Probably this FITS file does not conform OGIP standard. "
        "Only OGIP FITS files supported.","chkislc")
        return False
    if HDU.header['HDUCLAS1']!='LIGHTCURVE':
        printwarn("The 'HDUCLAS1' keyword is not 'LIGHTCURVE'. "
        "Probably this is not a light curve extension. ","chkislc")
        return False
    if not 'TIME' in HDU.columns.names:
        printwarn("There is no column 'TIME'.","chkislc")
        return False
    if not 'RATE' in HDU.columns.names:
        printwarn("There is no column 'RATE'.","chkislc")
        return False
    return True
    
def chkisgti(HDU):
    if not chkkey(HDU,'HDUCLASS'):
        printwarn("There is no 'HDUCLASS' keyword in the HDU header. "
        "Probably this FITS file does not conform OGIP standard. "
        "Only OGIP FITS files supported.","chkisgti")
        return False
    if HDU.header['HDUCLASS']!='OGIP':
        printwarn("'HDUCLASS' keyword is not 'OGIP'. "
        "Only OGIP FITS files supported.","chkisgti")
        return False
    if not chkkey(HDU,'HDUCLAS1'):
        printwarn("There is no 'HDUCLAS1' keyword in the HDU header. "
        "Probably this FITS file does not conform OGIP standard. "
        "Only OGIP FITS files supported.","chkisgti")
        return False
    if HDU.header['HDUCLAS1']!='GTI':
        printwarn("The 'HDUCLAS1' keyword is not 'GTI'. "
        "Probably this is not a GTI extension. ","chkisgti")
        return False
    if not 'START' in HDU.columns.names:
        printwarn("There is no column 'START'.","chkisgti")
        return False
    if not 'STOP' in HDU.columns.names:
        printwarn("There is no column 'STOP'.","chkisgti")
        return False
    return True
    
def name_generator(refname, mod, only_name=False):
  '''Generate new filenames using refrence name and modificators'''
  
  check_keywords_are_allowed(mod,['starts_with','ends_with','add_extension','replace_extension'])
  
  res = os.path.basename(refname)       
  if 'starts_with' in mod:
    res = mod['starts_with']+res
  if 'ends_with' in mod:
    res_split = os.path.splitext(res)
    res ='{}{}{}'.format(res_split[0],mod['ends_with'],res_split[1])  #Add '_xxx' before extension
  if 'add_extension' in mod:
    res = res + '.' + mod['add_extension']
  if 'replace_extension' in mod:
    res = os.path.splitext(res)[0]+'.' + mod['replace_extension']
    
  if only_name or not os.path.dirname(refname):
    return res
  else:
    res = os.path.dirname(refname) + '/' + res
  return res 
    
    
#############################  DS9  ####################################

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


##################### XSELECT ###############################
def __xselect_command_generator(data, header=None):
  '''Generate extractions commands. The 'data' is a list of dicts. Each dict contains
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
    data_for_generator[0]['extra_commands']=['clear region']    #Remove object's region before background extractions
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

  
###################################### XSPEC ###########################
def xspec_init(set_rebin=(5,5), statMethod='cstat', show_background=False, xaxis='keV', query='yes'):
  import xspec
  xspec.Plot.xAxis = xaxis
  xspec.Plot.background=show_background
  xspec.Plot.setRebin(*set_rebin)
  xspec.Fit.query=query
  xspec.Fit.statMethod = statMethod
  return xspec
  
def xspec_load_spectra(xspec, data_struct, clear=True, quiet=False):
  ''' data_struct must be a list of dicts '''
  
  
  #TODO make data_struct a class
  if type(data_struct) is not list:   
    raise TypeError("'data_struct' must be a list of dicts")
  
  #TODO: relative paths of backfile, respfile etc. may cause conflicts with change_dir()
  res=[]  #Return list of xspec.Spectrum() objects
  startWD=os.getcwd()
  startChatter=xspec.Xset.chatter
  
  if clear:
    xspec.AllData.clear()
    
  for entry in data_struct:
    xspec.Xset.chatter=0  #Hide the output because it shows message with no bkg and resp files which may be confusing
    check_keywords_are_allowed(entry,['specfile','bkgfile','respfile','arffile','ignore'])
    
    #Go the spectrum host dir to load correct bkg and resp files from its FITS-header
    os.chdir(os.path.dirname(entry['specfile'])) 
    sp = xspec.Spectrum(entry['specfile'])
    if 'bkgfile' in entry:
      sp.background=entry['bkgfile']
      sp.show()
    if 'respfile' in entry:
      sp.response=entry['respfile']
    if 'arffile' in entry:
      sp.response.arf=entry['arffile']
    if 'ignore' in entry:
      for line in entry['ignore']:
        sp.ignore(line)
    
    os.chdir(startWD)
    res.append(sp)
    xspec.Xset.chatter=startChatter
    if not quiet:
      sp.show()
    
  return res
  
def xspec_load_single_spectrum(xspec, specfile, bkgfile='', respfile='', arffile='', ignore=[], clear=True, quiet=False):
  '''Append new spectrum to already loaded'''

  data_struct = {'specfile':specfile}
  if bkgfile:
    data_struct['bkgfile']=bkgfile
  if respfile:
    data_struct['respfile']=respfile
  if arffile:
    data_struct['arffile']=arffile
  if ignore:
    data_struct['ignore']=ignore
    
  return xspec_load_spectra(xspec, [data_struct], clear, quiet)[0]

  
  
########################################################################
def group_spectrum_with_grppha(infile, outfile='', group_min_counts=1, clobber=False):
  
  if not outfile:
    outfile=name_generator(infile,{'ends_with':'_grp{}'.format(group_min_counts)})   #Add _grpN to the infile name
  
  check_file_not_exist_or_remove(outfile, overwrite=clobber,action='exception')  
  if callftools('grppha {} {}'.format(infile,outfile), stdin='group min {}\nexit\n'.format(group_min_counts)):
    return outfile

  return None

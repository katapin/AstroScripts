#!/usr/bin/python
"""Functions to perform some routine procedures involving ds9."""

import pyds9
import random
from numpy import ndarray
import os, sys, argparse



#### Defaults ####
ARG_DEFAULTS_CMAP='bb'
ARG_DEFAULTS_ZOOMTO='fit'
ARG_DEFAULTS_SCALE='linear'
ARG_DEFAULTS_SCALEMODE='zscale'
ARG_DEFAULTS_TITLE='NBI'


#### Conent from base.py

def ask_user(question:str, allowed_answers:list=[''], greetings:str=None):
    if greetings: print(greetings)
    while True:
        answ=input(question)
        if answ in allowed_answers:
            return answ

def check_file_exists_or_die(filepath:str):
    """Check the file exists and die if it doesn't, returns absolute path."""
    if not os.path.exists(filepath):
        die(f"File '{filepath}' is not found")
    return os.path.abspath(filepath)
    
def check_file_is_FITS(filepath:str):
    """Check the file exists and die if it doesn't, returns absolute path."""
    from astropy.io import fits
    try:
        with fits.open(filepath):
            pass
    except OSError:
        return False
    return True    

def read_ascii_filelist(ascii_path:str, comment_symbol:str='#') -> list:
    """Read ascii file ignoring comments."""
    files=[]
    with open(ascii_path,'r') as fl:    #Read lines from file
        for line in fl:
            stripped = line.strip()
            if not stripped.startswith(comment_symbol):  #Remove lines starting with #
                word = stripped.split(comment_symbol)[0]
                #Now we'll try to check that the file is really a file list and
                # not a random file with text
                if len(word.split(' '))>1:
                    die("Wrong format of ASCII filelist. Filenames mustn't contain "\
                        f"spaces but string '{word}' does.")
                files.append(check_file_exists_or_die(word.strip())) #Remove comments at the end of the line 
    return files


__COLOR_SEQUENCES={'red':'91m','yellow':'93m','green':'92m',
'cyan':'96m', 'bold':'01m'}

def _print_colored(text, color):
    print('\033['+__COLOR_SEQUENCES[color]+text+'\033[0m')
    
def printerr(text:str):
    """Print message in red."""
    _print_colored('Error: '+text, 'red')
    
def printwarn(text:str):
    """Print message in yellow."""
    _print_colored('Warning: '+text, 'yellow')
    
def printbold(text:str):
    """Print message with bold font."""
    _print_colored(text, 'bold')


def die(text:str, code:int=1):
    """Print error message and exit."""
    printerr(text)
    sys.exit(code)
    

####  Exceptions for ds9 ##

class DS9Error(Exception):
    """Exception to handle errors caused by ds9 program."""
    def __init__(self, ds9obj:pyds9.DS9, command:str, custom_message:str=None):
        self.ds9obj = ds9obj
        self.command=command
        if custom_message is not None:
            self.msg = custom_message
        else:
            self.msg = f"An error occured in ds9 title='{ds9obj.target}' during "\
                f"processing the command '{command}'."
        super().__init__(self.msg)

class DS9CommandUndefiedError(DS9Error):
    """Exeption to raise when the ds9 can't recognize the XPA command send."""
    def __init__(self, ds9obj:pyds9.DS9, command:str):
        msg = f"The ds9 title='{ds9obj.target}' can't recognize "\
            f"command '{command}'."
        super().__init__(ds9obj, command, msg)
                
class DS9WindowClosedError(DS9Error):
    """Exeption to raise when the ds9 process turned out to be closed."""
    def __init__(self, ds9obj:pyds9.DS9, command:str):
        msg=f"The ds9 window title='{ds9obj.target}' was closed."
        super().__init__(ds9obj, command, msg)
        
class DS9FileLoadError(DS9Error):
    """Exeption to raise when a file can't be loaded."""
    def __init__(self, ds9obj:pyds9.DS9, command:str, filepath:str):
        msg=f"The ds9 title='{ds9obj.target}' can't load file '{filepath}'"
        self.filepath = filepath
        super().__init__(ds9obj, command, msg)
    

        
#### ds9 functions

def start(filename:str=None, *, workdir:str = None, title:str = None, wait:int=20, 
             cmap:str=None, scale:str=None, scalemode:str=None, zoomto:str = None, 
             extra_commands:list = None):
    
    if title is None:       #Generate random title for new ds9 window
        title='ds9_'+'{:04d}'.format(random.randint(0,9999))
    ds9obj = pyds9.DS9(target=title, wait=wait)
      
    commands=[]
    
    if workdir: send_command(ds9obj, f'cd {workdir}')
    if filename: send_command(ds9obj, f'file {filename}')
    if zoomto: send_command(ds9obj, f'zoom to {zoomto}')
    if scale: send_command(ds9obj, f'scale {scale}')
    if scalemode: send_command(ds9obj, f'scale mode {scalemode}')
    if cmap: send_command(ds9obj, f'cmap {cmap}')
    
    for cmd in (commands + (extra_commands or [])): 
        send_command(ds9obj, cmd)
    return ds9obj
  
def load_region(ds9obj, regfile, color:str=None, group:str=None, clear:bool=False):

    if clear:
         send_command(ds9obj, 'region delete')
      
    buf=None  #Data buffer for pyds9.set command    
    if color or group:  #This demands to edit the region file content
        with open(regfile,'r') as f:
            content=f.readlines()
        extra_text=' #'
        if color:
            extra_text+=' color={}'.format(color)
        if group:
            extra_text+=' group={}'.format(color)
          
        for i in range(3,len(content)):   #Ignore three first header lines 
            content[i]=content[i].replace('\n',extra_text+'\n')   #Add to each string of 
            
        cmd='regions'
        buf=''.join(content)
    else:
        cmd='region load {}'.format(regfile)
        
    send_command(ds9obj, cmd, buf)
            

def send_command(ds9obj:pyds9.DS9, command:str, buf=None):
    """Send command to the ds9 process and raise exception if error occures."""
    try:
        # print(command)
        if ds9obj.set(command, buf=buf) == 0:
            raise DS9Error(ds9obj, command)
    except ValueError as ex:
        if ex.args[0].startswith('ds9 is no longer running'):
            raise DS9WindowClosedError(ds9obj, command) from ex
        elif ex.args[0].startswith('XPA$ERROR Unable to load'):
            filename = ex.args[0].split(' ')[6]
            raise DS9FileLoadError(ds9obj, command, filename) from ex
        elif ex.args[0].startswith('XPA$ERROR undefined command'):
            raise DS9CommandUndefiedError(ds9obj, command) from ex
        else:
            raise DS9Error(ds9obj, command) from ex

  
def load_image(ds9obj, image, slot:int=None, *, regfile:str=None, zoomto:str=None,
       extra_commands:list=None):
    
    if slot is not None:
        frame_list = ds9obj.get('frame all').split(' ')  #Make list from string  '1 2 3'
        if str(slot) not in frame_list:
            raise DS9Error(f"Slot '{slot}' is not found", ds9obj=ds9obj)
        send_command(ds9obj, f'frame {slot}')
        
    if isinstance(image, str):  #it's a file path
        send_command(ds9obj, f'file {image}')
    elif isinstance(image, ndarray):
        raise NotImplementedError()

    if regfile: load_region(ds9obj, regfile)
    if zoomto: send_command(ds9obj, f'zoom to {zoomto}')
    if extra_commands: 
        for cmd in extra_commands: send_command(ds9obj, cmd)
        
def show_iteratively(ds9obj, filelist:list, *, preprocessor=None, asker=None,
        region:str=None, zoomto:str=None, extra_per_file:str=None):
    n=len(filelist)
    for i,file in enumerate(filelist):
        print(f'#{i+1}/{n} {file}:')
        try:
            load_image(ds9obj, file, regfile=region, zoomto=zoomto,
                       extra_commands=extra_per_file)
            asker(file)
        except DS9FileLoadError as ex:
            printwarn(f"Can't load '{ex.filepath}'. May be it's not an image.")
    
    
def _main():
    parser = argparse.ArgumentParser(description="Iteratively shows FITS images in ds9")
    parser.add_argument('input', nargs='+', help="One or multiple FITS images or ASCII list")
    parser.add_argument('--cmap', nargs='?', help=f"Name of the colormap in ds9. Default is '{ARG_DEFAULTS_CMAP}'.",
        default=ARG_DEFAULTS_CMAP)
    parser.add_argument('--scale', nargs='?', help=f"Scale parameter of ds9. Default is '{ARG_DEFAULTS_SCALE}'.",
        default=ARG_DEFAULTS_SCALE)
    parser.add_argument('--scale-mode', nargs='?', help=f"Scale mode. Default is '{ARG_DEFAULTS_SCALEMODE}'.",
        default=ARG_DEFAULTS_SCALEMODE)
    parser.add_argument('--zoom-to', nargs='?', help=f"Scale parameter of ds9. Default is '{ARG_DEFAULTS_ZOOMTO}'.",
        default=ARG_DEFAULTS_ZOOMTO)
    parser.add_argument('--region', nargs='?', help="Load region file for every image")
    parser.add_argument('--title', nargs='?', help="Title of the ds9 window. It allows "
        "to connect to the existing session. The title will be random if "
        f"this argument is empty. Default is '{ARG_DEFAULTS_TITLE}'.", default=ARG_DEFAULTS_TITLE)
    parser.add_argument('--ds9-once', nargs='*', help="Extra arguments for ds9, run at start "\
        "(before any file will be loaded)", action='append', default=[])
    parser.add_argument('--ds9-per-file', nargs='*', help="Extra arguments for ds9, run for "\
        "each file (after the file will be loaded).", action='append', default=[])
    
    argnspace=parser.parse_args(sys.argv[1:])
    
    if len(argnspace.input)==1: #try to guess type of the only infile 
        infile = check_file_exists_or_die(argnspace.input[0])
        if check_file_is_FITS(infile):  #The file is a FITS
            filelist=[infile]
        else:                                    #The file is an ASCII list
            try: 
                filelist=read_ascii_filelist(infile)
            except UnicodeError:  #It's not an ASCII list...
                die(f"Can't determine type of the input file '{infile}'")
            if len(filelist) == 0:
                print(f"The list '{infile}' is empty. Noting to do. Exiting...")
                sys.exit(0)
    else:
        filelist=argnspace.input
    
    #Run ds9 process (or connect to existing)
    commands_once = [ item for sublist in argnspace.ds9_once for item in sublist ] 
    ds9obj=start(title=argnspace.title, workdir=os.getcwd(), cmap=argnspace.cmap, 
             scale=argnspace.scale, scalemode=argnspace.scale_mode, extra_commands=commands_once)
    
    
    commands_perfile = [ item for sublist in argnspace.ds9_per_file for item in sublist ]
    def asker(x): ask_user("Press 'Enter' to continue")
    show_iteratively(ds9obj, filelist, asker=asker, region=argnspace.region, 
             zoomto=argnspace.zoom_to, extra_per_file=commands_perfile)

if __name__ == '__main__':    
    try:
        _main()
    except Exception as ex:
        die(ex.args[0])
    # _main()
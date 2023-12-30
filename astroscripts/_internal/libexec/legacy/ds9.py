#!/usr/bin/python
"""Functions to perform some routine procedures involving ds9."""

import pyds9
import random
from numpy import ndarray

class DS9Error(Exception):
    """Exception to handle errors caused by ds9 program."""
    def __init__(self, ds9obj:pyds9.DS9, custom_message:str=None):
        self.msg = f"An error occured in ds9 title='{ds9obj.target}'."
        super().__init__(self.msg)

class DS9CommandUndefiedError(DS9Error):
    
    def __init__(self, ds9obj:pyds9.DS9, cmd:str, custom_message:str=None):
        self.ds9obj = ds9obj
        self.cmd=cmd
        if custom_message:
            self.msg = custom_message
        else:
            self.msg = "An error occured in ds9 title='{}' during "\
                "processing command '{}'.".format(ds9obj.target, cmd)
        super().__init__(self.msg)
        
                
class DS9WindowClosedError(DS9Error):
    """Exeption to raise when the ds9 process turned out to be closed."""
    def __init__(self, ds9obj:pyds9.DS9, cmd:str, custom_message:str=None):
        msg=f"The ds9 window title='{ds9obj.target}' was closed."
        super().__init__(ds9obj, cmd, msg)
        
class DS9FileLoadError(DS9Error):
    pass
    


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
            

def send_command(ds9obj:pyds9.DS9, cmd:str, buf=None):
    """Send command to the ds9 process and raise exception if error occures."""
    try:
        print(cmd)
        if ds9obj.set(cmd, buf=buf) == 0:
            raise DS9Error(ds9obj, cmd)
    except ValueError as ex:
        print(ex)
        if ex.args[0].startswith('ds9 is no longer running'):
            raise DS9Closed(ds9obj, cmd) from ex
        elif
            'XPA$ERROR Unable to load'
        elif 'XPA$ERROR undefined command'        
    
            
        raise DS9Error(ds9obj, cmd) from ex

  
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

    if regfile is not None: load_region(ds9obj, regfile)
    if zoomto is not None: send_command(ds9obj, f'zoom to {zoomto}')
    if extra_commands: 
        for cmd in extra_commands: send_command(ds9obj, cmd)
    
    
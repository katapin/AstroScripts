#!/usr/bin/python

import sys
import os
import xml.etree.ElementTree as ET
import lxml.etree as lxmlet
import re


def die(s):
    print("Error: " + s)
    sys.exit(1);



class Dataset:
	def __init__(self, expid, duration, expmode):
		self.expid=expid
		self.duration=duration
		self.mode=expmode
		
		##Check mode: Is it imaging mode? 
		#Imaging modes are 'Prime...'	skip=0
		#Timing modes are 'Fast...'		slip=1
		#Full list of modes
		#http://xmm2.esac.esa.int/external/xmm_sw_cal/calib/documentation/CALHB/node769.html
	
		self.skip=0 if (re.match('^Prime*',expmode)) else 1
	
class Observation:
	def __init__(self):
	#
	##Init class members
	
		self.obsid=""
		self.source=""
		self.ra=""
		self.dec=""
		self.expnum={"pn":-1, "mos1":-1, "mos2":-1}
		self.data={"pn":[], "mos1":[], "mos2":[]}
			
		
	def __parseexposures(self, INSTR, key):
	# 
	##Explore every exposure, check mode, duration and id.
	
		explist=[]
		xmlexplist=INSTR.findall("EXPOSURE")
		self.expnum[key]=len(xmlexplist) 	#True number of exposures
		for EXP in xmlexplist:
			expmode=EXP.get("mode")
			expid=EXP.get("expid")
			duration=EXP.get("duration")
			curdataset=Dataset(expid,duration,expmode)
			explist.append(curdataset)
			self.data[key].append(curdataset)

		
	def loadfromxmmsas(self,xmlfilename):
	#
	##Load xmmextractor config and parse with typical structure:

#<BODY> 
#    <CONFIG>
#        <OBSERVATION id=... dafault=...>
#        ...
#        </OBSERVATION>
#        <INSTRUMENT value=...>
#            <EXPOSURE expid=... duration=...>
#            ...
#            </EXPOSURE>     
#        </INSTRUMENT
	##
	#
		
		#Check input
		if not os.path.isfile(xmlfilename):
			die("XML-file is not found.")
		
		#Load xml
		tree = ET.parse(xmlfilename)
		root = tree.getroot() 
		if (root.tag!="BODY"): 
			die("Can't find tag 'BODY'. Wrong xml format. ")
		if (root[0][0].tag!="OBSERVATION"):
			die("Can't find tag 'OBSERVATION'. Wrong xml format. ")

		#Load info 
		OBS=root[0][0]		#OBSERVATION
		for PARAM in OBS:
			parname=PARAM.get("id")
			if parname=="obsid":
				self.obsid=PARAM.get("default")
			elif parname=="sourcename":
				self.source=PARAM.get("default")
			elif parname=="ra":
				self.ra=PARAM.get("default")
			elif parname=="dec":
				self.dec=PARAM.get("default")
			elif parname=="EPN":
				self.expnum["pn"] = 0 if PARAM.get("default")=="no" else 1
			elif parname=="EMOS1":
				self.expnum["mos1"] = 0 if PARAM.get("default")=="no" else 1
			elif parname=="EMOS2":
				self.expnum["mos2"] = 0 if PARAM.get("default")=="no" else 1

		#Check acquired information
		if (self.obsid==""):
			die("Can't find ObsID number. Wrong xml format.")
		if (self.expnum["pn"]==-1):
			die("Can't find information about EPIC-PN. Wrong xml format.")
		if (self.expnum["mos1"]==-1):
			die("Can't find information about EPIC-MOS1. Wrong xml format.")
		if (self.expnum["mos2"]==-1):
			die("Can't find information about EPIC-MOS2. Wrong xml format.")
			
		
		for INSTR in root[0].findall("INSTRUMENT"):
			if INSTR.get("value")=="EPN":
				if self.expnum["pn"]: self.__parseexposures(INSTR,"pn")
			elif INSTR.get("value")=="EMOS1":
				if self.expnum["mos1"]: self.__parseexposures(INSTR,"mos1")
			elif INSTR.get("value")=="EMOS2":
				if self.expnum["mos2"]: self.__parseexposures(INSTR,"mos2")
	
	def savexml(self,filename):
	#
	##Save class date into xml-file
		
		obs=ET.Element('OBSERVATION')
		ET.SubElement(obs,'PARAM',name="obsid",value=self.obsid)
		ET.SubElement(obs,'PARAM',name="sourcename",value=self.source)
		ET.SubElement(obs,'PARAM',name="ra",value=self.ra)
		ET.SubElement(obs,'PARAM',name="dec",value=self.dec)
		ET.SubElement(obs,'PARAM',name="expnum_pn",value=str(self.expnum['pn']))
		ET.SubElement(obs,'PARAM',name="expnum_mos1",value=str(self.expnum['mos1']))
		ET.SubElement(obs,'PARAM',name="expnum_mos2",value=str(self.expnum['mos2']))
		
		pndata=ET.SubElement(obs,'PNDATA')
		mos2data=ET.SubElement(obs,'MOS1DATA')
		mos1data=ET.SubElement(obs,'MOS2DATA')
		instr_branch={'pn':pndata, 'mos1':mos1data, 'mos2':mos2data}
		instr_list=['pn', 'mos1', 'mos2']
		for j in range(0,3):
			key=instr_list[j]
			for i in range(0,self.expnum[key]):
				curdataset=ET.SubElement(instr_branch[key],'DATASET')
				ET.SubElement(curdataset,'PARAM',name="expid",value=str(self.data[key][i].expid))
				ET.SubElement(curdataset,'PARAM',name="mode",value=str(self.data[key][i].mode))
				ET.SubElement(curdataset,'PARAM',name="duration",value=str(self.data[key][i].duration))
		
		
		
		tree=ET.ElementTree(obs)
		tree.write("my.xml", encoding='utf-8', xml_declaration=True) 
		
		print(ET.tostring(obs))
		#osb.
		#tree=ET.ElementTree()
		#root=tree.getroot()
		##root.append('a')
		#ET.SubElement(root,'a')
		#tree.SubElement('a')
		#print(tree)
		
		#tree=lxmlet.Element("OBSERVATION")
		#lxmlet.SubElement(tree, "a")
		#print(lxmlet.tostring(tree))
		#tree=ET.fromstring("<?xml version='1.0'><QWE></QWE>")
		##tree=ET()
		
		#tree=ET.parse('0144rundefault.xml')
	#	tree.write('my.xml')


##Chech initial data
#PWD=os.getcwd()
#if not os.environ.get('SAS_ODF'):
	#die("Variable 'SAS_ODF' is not defined. Please define the variable and try again.")
#else:
	#ODF=os.environ.get('SAS_ODF')
	#if not os.path.isfile(ODF):
		#die("Variable 'SAS_ODF does not point to the valid summary-file.")

#if not os.environ.get('SAS_CCF'):
	#die("Variable 'SAS_CCF' is not defined. Please define the variable and try again.")
#else:
	#CCF=os.environ.get('SAS_CCF')
	#if not os.path.isfile(CCF):
		#die("Variable 'SAS_CCF does not point to the valid CCF-file.")
   

			
			
obs = Observation()
obs.loadfromxmmsas('0144rundefault.xml')
obs.savexml('my.xml')


#print(root.tag)
#print(OBS)
#print(root[0][1].tag)


#! /usr/bin/env python


import os
import sys
import uuid
from math import pow
from string import maketrans

def enum(*sequential, **named):
   enums = dict(zip(sequential, range(len(sequential))), **named)
   reverse = dict((value, key) for key, value in enums.iteritems())
   enums['reverse_mapping'] = reverse
   return type('Enum', (), enums)

class UiElement:
   name = ''
   ident = ''
   description = None
   parent = None
   children = None


class Group(UiElement):
   collapse = True
   

   def __init__(self, name, collapse=True, description=None, ident=None):
      self.name = name
      self.collapse = collapse
      self.description = description
      
      if ident==None:
         ident = name.replace(' ', '')
         ident = ident.replace('_', '')
         ident = ident.lower()
         temptrans = maketrans('','')
         self.ident = ident.translate(temptrans,'aeiouy')
         
      else:
         self.ident = ident

      self.unique_name = None

   def __str__(self):
      return 'Group %s' % self.name


class Tab(Group):
   def __init__(self, name, collapse=True, description=None,  ident=None):
      self.name = name
      self.collapse = collapse
      self.description = description
      
      if ident==None:
         ident = name.replace(' ', '')
         ident = ident.replace('_', '')
         ident = ident.lower()
         temptrans = maketrans('','')
         self.ident = ident.translate(temptrans,'aeiouy')
      else:
         self.ident = ident

      self.unique_name = None

   def __str__(self):
      return 'Tab %s' % self.name

class Parameter(UiElement):
   ptype = None
   label = None
   mn = None
   mx = None
   smn = None
   smx = None
   connectible = True
   enum_names = None
   default = None
   fig = None
   figc = None
   short_description = ''
   presets = None

   def __init__(self, name, ptype, default, label=None, description=None, mn=None, mx=None, smn=None, smx=None, connectible=True, enum_names=None, fig=None, figc=None, presets={}):
      self.name = name
      self.ptype = ptype
      self.default = default
      self.description = description
      self.presets = presets

      if description is not None:
         if len(description) > 150:
            # build the short description, < 200 chars
            sentences = description.split('.')
            sdlen = 0
            for s in sentences:
               self.short_description += s + "."
               sdlen += len(s)
               if sdlen > 150:
                  break
         else:
            self.short_description = description

      if not label:
         label = name   
      self.label = label   
      self.mn = mn
      self.mx = mx
      self.smn = smn
      self.smx = smx
      self.fig = fig
      self.figc = figc

      if ptype not in ('bool', 'int', 'float', 'rgb', 'vector', 'string', 'enum', 'matrix'):
         raise ValueError('parameter %s has unrecognized ptype %r' % (name, ptype))

      if ptype == 'enum':
         self.enum_names = enum_names

      if ptype == 'bool' or ptype == 'string' or ptype == 'int' or ptype == 'enum':
         self.connectible = False
      else:
         self.connectible = connectible

      # do a quick sanity check on the ptype and default
      if ptype == 'bool' and type(default) is not bool:
         raise ValueError('parameter %s was typed as %s but default had type %r' % (name, ptype, type(default)))
      elif ptype in ('rgb', 'vector') and type(default) is not tuple:
         raise ValueError('parameter %s was typed as %s but default had type %r' % (name, ptype, type(default)))
      elif ptype in ('string', 'enum') and type(default) is not str:
         raise ValueError('parameter %s was typed as %s but default had type %r' % (name, ptype, type(default)))

      self.default = default

   def __str__(self):
      return 'Parameter %s %s' % (self.name, self.ptype)

class AOV(Parameter):
   def __init__(self, name, ptype, label=None):
      self.name = name
      self.ptype = ptype
      self.label = label
      self.connectible = False

   def __str__(self):
      return 'AOV %s %s' % (self.name, self.ptype)

class ShaderDef:
   name = None
   intro = None
   description = None
   output = None
   c4d_classification = None
   c4d_menu = None
   c4d_command_id = None
   maya_name = None
   maya_classification = None
   maya_id = None
   maya_swatch = False
   maya_bump = False
   maya_matte = False
   soft_name = None
   soft_classification = None
   help_url = None
   houdini_icon = None
   soft_version = 1
   hierarchy_depth = 0
   root = None
   current_parent = None
   parameters = []
   aovs = []
   groups = []
   tabs = []

   def __init__(self):
      self.root = Group('__ROOT__', False)
      self.current_parent = self.root

   def setMayaClassification(self, mc):
      self.maya_classification = mc

   def setDescription(self, d):
      self.description = d

   def setOutput(self, o):
      self.output = o

   def beginGroup(self, name, collapse=True, description='', ident=None):
      g = None
      if self.current_parent is self.root:
         g = Tab(name, collapse, description,ident)
      else:
         g = Group(name, collapse, description,ident)

      if not self.current_parent.children:
         self.current_parent.children = [g]
      else:
         self.current_parent.children.append(g)

      g.parent = self.current_parent
      self.current_parent = g

   def endGroup(self):
      self.current_parent = self.current_parent.parent

   def parameter(self, name, ptype, default, label=None, description=None, mn=None, mx=None, smn=None, smx=None, connectible=True, enum_names=None, fig=None, figc = None, presets=None):
      p = Parameter(name, ptype, default, label, description, mn, mx, smn, smx, connectible, enum_names, fig, figc, presets)
      if not self.current_parent.children:
         self.current_parent.children = [p]
      else:
         self.current_parent.children.append(p)

      p.parent = self.current_parent
      self.parameters.append(p)

   def aov(self, name, ptype, label=None, description=None):
      p = AOV(name, ptype, label)
      if not self.current_parent.children:
         self.current_parent.children = [p]
      else:
         self.current_parent.children.append(p)

      p.parent = self.current_parent
      self.aovs.append(p)

   def __str__(self):
      return 'ShaderDef %s' % self.name

   def shader(self, d):
      self.name = d['name']
      self.intro = d['intro']
      self.description = d['description']
      self.output = d['output']
      self.c4d_classification = d.get('c4d_classification')
      self.c4d_menu = d.get('c4d_menu')
      self.c4d_command_id = d.get('c4d_command_id')
      self.maya_name = d['maya_name']
      self.maya_classification = d['maya_classification']
      self.maya_id = d['maya_id']
      self.maya_swatch = d['maya_swatch']
      self.maya_bump = d['maya_bump']
      self.maya_matte = d['maya_matte']
      self.soft_name = d['soft_name']
      self.soft_classification = d['soft_classification']
      self.soft_version = d['soft_version']
      if 'help_url' in d.keys():
         self.help_url = d['help_url']
      else:
         self.help_url = 'https://bitbucket.org/anderslanglands/alshaders/wiki/Home'
      if 'houdini_icon' in d.keys():
         self.houdini_icon = d['houdini_icon']


   #all groups need unique ident in houdini
   def uniqueGroupIdents(self):
      grpidents = []    
      for grp in self.groups:
         oldident=grp.ident
         if grp.ident in grpidents:
            grp.ident = '%s%d' % (grp.ident, grpidents.count(oldident))
         grpidents.append(oldident)    

class group:
   name = ''
   collapse = True
   description = ''
   ui = None
   ident = None

   def __init__(self, ui, name, collapse=True, description=None, ident=None):
      self.ui = ui
      self.name = name     
      self.collapse = collapse
      self.description = description
      self.ident = ident

   def __enter__(self):
      self.ui.beginGroup(self.name, self.collapse, self.description, self.ident)

   def __exit__(self, type, value, traceback):
      self.ui.endGroup()


def DebugPrintElement(el, d):
   indent = ''
   for i in range(d):
      indent += '\t'
   print '%s %s' % (indent, el)
   if isinstance(el, Group) and el.children:
      print '%s has %d children' % (el, len(el.children))
      for e in el.children:
         DebugPrintElement(e, d+1)
      d -= 1

def DebugPrintShaderDef(sd):
   print '%s' % sd
   DebugPrintElement(sd.root, 0)

def writei(f, s, d=0):
   indent = ''
   for i in range(d):
      indent += '\t'
   f.write('%s%s\n' %(indent, s))

def WalkAETemplate(el, f, d):
   if isinstance(el, Group):
      writei(f, 'self.beginLayout("%s", collapse=%r)' % (el.name, el.collapse), d)

      if el.children:
         for e in el.children:
            WalkAETemplate(e, f, d)

      writei(f, 'self.endLayout() # END %s' % el.name, d)

   elif isinstance(el, Parameter):
      if el.ptype == 'float':
         writei(f, 'self.addCustomFlt("%s")' % el.name, d)
      elif el.ptype == 'rgb':
         writei(f, 'self.addCustomRgb("%s")' % el.name, d)
      else:
         writei(f, 'self.addControl("%s", label="%s", annotation="%s")' % (el.name, el.label, el.short_description), d)

def WalkAETemplateParams(el, f, d):
   if isinstance(el, Group):
      if el.children:
         for e in el.children:
            WalkAETemplateParams(e, f, d)

   elif isinstance(el, Parameter):
      writei(f, 'self.params["%s"] = Param("%s", "%s", "%s", "%s", presets=%s)' % (el.name, el.name, el.label, el.short_description, el.ptype, el.presets), d)

def WriteAETemplate(sd, fn):
   f = open(fn, 'w')
   writei(f, 'import pymel.core as pm', 0)
   writei(f, 'from alShaders import *\n', 0)

   writei(f, 'class AE%sTemplate(alShadersTemplate):' % sd.name, 0)
   writei(f, 'controls = {}', 1)
   writei(f, 'params = {}', 1)
   
   writei(f, 'def setup(self):', 1)

   writei(f, 'self.params.clear()', 2)
   for e in sd.root.children:
      WalkAETemplateParams(e, f, 2)
   
   writei(f, '')

   if sd.maya_swatch:
      writei(f, 'self.addSwatch()', 2)

   writei(f, 'self.beginScrollLayout()', 2) # begin main scrollLayout
   writei(f, '')

   if sd.maya_matte:
      writei(f, 'self.beginLayout("Matte")', 2)
      writei(f, 'self.addControl("aiEnableMatte", label="Enable Matte")', 2)
      writei(f, 'self.addControl("aiMatteColor", label="Matte Color")', 2)
      writei(f, 'self.addControl("aiMatteColorA", label="Matte Opacity")', 2)
      writei(f, 'self.endLayout()', 2)


   for e in sd.root.children:
      WalkAETemplate(e, f, 2)

   if sd.maya_bump:
      writei(f, 'self.addBumpLayout()', 2)      

   writei(f, '')
   writei(f, 'pm.mel.AEdependNodeTemplate(self.nodeName)', 2)

   writei(f, 'self.addExtraControls()', 2)

   writei(f, '')
   writei(f, 'self.endScrollLayout()', 2) #end main scrollLayout

   f.close()

#########
# C4DtoA
#########

# Parameter ids are generated from the name.
def GenerateC4DtoAId(name, parameter_name):
   unique_name = "%s.%s" % (name, parameter_name)

   pid = 5381
   for c in unique_name:
      pid = ((pid << 5) + pid) + ord(c) # hash*33 + c
      #pid = pid*33 + ord(c)
      #print "%d %c" % (pid, c)

   # convert to unsigned int
   pid = pid & 0xffffffff
   if pid > 2147483647: pid = 2*2147483647 - pid + 2

   return pid

# Writes out group ids to the header file.
def WalkC4DtoAHeaderGroups(el, f, name):
   if isinstance(el, Group):
      group_name = el.unique_name if el.unique_name else el.name
      writei(f, 'C4DAI_%s_%s_GRP,' % (name.upper(), group_name.upper().replace(' ', '_')), 1)

      if el.children:
         for e in el.children:
            WalkC4DtoAHeaderGroups(e, f, name)

# Writes out parameter ids to the header file.
def WalkC4DtoAHeaderParameters(el, f, name):
   if isinstance(el, Group):
      if el.children:
         for e in el.children:
            WalkC4DtoAHeaderParameters(e, f, name)

   # skip aovs
   #if isinstance(el, AOV):
   #  pass

   if isinstance(el, Parameter):
      pid = GenerateC4DtoAId(name, el.name)
      writei(f, 'C4DAIP_%s_%s = %d,' % (name.upper(), el.name.upper().replace(' ', '_'), pid), 1)

# Writes .h header file.
def WriteC4DtoAHeaderFile(sd, name, build_dir):
   path = os.path.join(build_dir, "C4DtoA", "res", "description", "ainode_%s.h" % name)
   if os.path.exists(path):
      os.remove(path)
   if not os.path.exists(os.path.dirname(path)):
      os.makedirs(os.path.dirname(path))
   f = open(path, 'w')

   writei(f, '#ifndef _ainode_%s_h_' % name, 0)
   writei(f, '#define _ainode_%s_h_' % name, 0) 
   writei(f, '', 0) 
   writei(f, 'enum', 0) 
   writei(f, '{', 0) 
   writei(f, 'C4DAI_%s_MAIN_GRP = 2001,' % (name.upper()), 1)

   # groups
   for e in sd.root.children:
      WalkC4DtoAHeaderGroups(e, f, name)
   # matte
   if sd.maya_matte:
      writei(f, 'C4DAI_%s_MATTE_GRP,' % name.upper(), 1)   

   writei(f, '', 0)

   # parameters
   for e in sd.root.children:
      WalkC4DtoAHeaderParameters(e, f, name)
   # matte
   if sd.maya_matte:
      writei(f, 'C4DAIP_%s_AIENABLEMATTE = %d,' % (name.upper(), GenerateC4DtoAId(name, 'aiEnableMatte')), 1)
      writei(f, 'C4DAIP_%s_AIMATTECOLOR = %d,' % (name.upper(), GenerateC4DtoAId(name, 'aiMatteColor')), 1)
      writei(f, 'C4DAIP_%s_AIMATTECOLORA = %d,' % (name.upper(), GenerateC4DtoAId(name, 'aiMatteColorA')), 1)

   writei(f, '};', 0) 
   writei(f, '', 0) 
   writei(f, '#endif', 0) 
   writei(f, '', 0) 

# Writes out layout to the resource file.
def WalkC4DtoARes(el, f, name, d):
   if isinstance(el, Group):
      group_name = el.unique_name if el.unique_name else el.name      
      writei(f, 'GROUP C4DAI_%s_%s_GRP' % (name.upper(), group_name.upper().replace(' ', '_')), d)
      writei(f, '{', d)

      if not el.collapse:
         writei(f, 'DEFAULT 1;', d+1)
         writei(f, '', 0)

      if el.children:
         for e in el.children:
            WalkC4DtoARes(e, f, name, d+1)

      writei(f, '}', d)
      writei(f, '', 0)

   # skip aovs
   #elif isinstance(el, AOV):
   #  pass

   elif isinstance(el, Parameter):
      writei(f, 'AIPARAM C4DAIP_%s_%s {}' % (name.upper(), el.name.upper().replace(' ', '_')), d)

# Writes .res resource file.
def WriteC4DtoAResFile(sd, name, build_dir):
   path = os.path.join(build_dir, "C4DtoA", "res", "description", "ainode_%s.res" % name)
   if os.path.exists(path):
      os.remove(path)
   if not os.path.exists(os.path.dirname(path)):
      os.makedirs(os.path.dirname(path))
   f = open(path, 'w')

   writei(f, 'CONTAINER AINODE_%s' % name.upper(), 0)
   writei(f, '{', 0)
   writei(f, 'NAME ainode_%s;' % name, 1)
   writei(f, '', 0)
   writei(f, 'INCLUDE GVbase;', 1)
   writei(f, '', 0)
   writei(f, 'GROUP C4DAI_%s_MAIN_GRP' % name.upper(), 1)
   writei(f, '{', 1)
   writei(f, 'DEFAULT 1;', 2)
   writei(f, '', 0)

   # matte
   if sd.maya_matte:
      writei(f, 'GROUP C4DAI_%s_MATTE_GRP' % name.upper(), 2)
      writei(f, '{', 2)
      writei(f, 'AIPARAM C4DAIP_%s_AIENABLEMATTE {}' % name.upper(), 3)
      writei(f, 'AIPARAM C4DAIP_%s_AIMATTECOLOR {}' % name.upper(), 3)
      writei(f, 'AIPARAM C4DAIP_%s_AIMATTECOLORA {}' % name.upper(), 3)
      writei(f, '}', 2)
      writei(f, '', 0)

   # groups and parameters
   for e in sd.root.children:
      WalkC4DtoARes(e, f, name, 2)

   writei(f, '}', 1)
   writei(f, '}', 0)
   writei(f, '') 

# Writes out group labels to the string file.
def WalkC4DtoAStringGroups(el, f, name):
   if isinstance(el, Group):
      group_name = el.unique_name if el.unique_name else el.name
      writei(f, 'C4DAI_%s_%s_GRP   "%s";' % (name.upper(), group_name.upper().replace(' ', '_'), el.name), 1)

      if el.children:
         for e in el.children:
            WalkC4DtoAStringGroups(e, f, name)

# Writes out parameter labels to the string file.
def WalkC4DtoAStringParameters(el, f, name):
   if isinstance(el, Group):
      if el.children:
         for e in el.children:
            WalkC4DtoAStringParameters(e, f, name)

   if isinstance(el, Parameter) and el.label:
      writei(f, 'C4DAIP_%s_%s   "%s";' % (name.upper(), el.name.upper(), el.label), 1)

# Writes .str string file.
def WriteC4DtoAStringFile(sd, name, build_dir):
   path = os.path.join(build_dir, "C4DtoA", "res", "strings_us", "description", "ainode_%s.str" % name)
   if os.path.exists(path):
      os.remove(path)
   if not os.path.exists(os.path.dirname(path)):
      os.makedirs(os.path.dirname(path))
   f = open(path, 'w')

   writei(f, 'STRINGTABLE ainode_%s' % name, 0)
   writei(f, '{', 0)
   writei(f, 'ainode_%s   "Arnold %s node";' % (name, name), 1) 
   writei(f, '', 0)
   writei(f, 'C4DAI_%s_MAIN_GRP   "Main";' % (name.upper()), 1) 

   # groups
   for e in sd.root.children:
      WalkC4DtoAStringGroups(e, f, name)
   # matte
   if sd.maya_matte:
      writei(f, 'C4DAI_%s_MATTE_GRP   "Matte";' % name.upper(), 1)

   writei(f, '', 0)

   # parameters
   for e in sd.root.children:
      WalkC4DtoAStringParameters(e, f, name)
   # matte
   if sd.maya_matte:
      writei(f, 'C4DAIP_%s_AIENABLEMATTE   "Enable matte";' % name.upper(), 1)
      writei(f, 'C4DAIP_%s_AIMATTECOLOR   "Matte color";' % name.upper(), 1)
      writei(f, 'C4DAIP_%s_AIMATTECOLORA   "Matte opacity";' % name.upper(), 1)

   writei(f, '}', 0)
   writei(f, '', 0)

# Writes resource files for C4DtoA.
# To describe the UI we need three resource files:
#  - .h header file which contains parameter ids
#  - .res resource file which contains widget layout
#  - .str string file which contains labels
def WriteC4DtoAResourceFiles(sd, name, build_dir):
   # NOTE group list is build in WriteMTD
   # create unique group names
   names = set()
   groups = sd.tabs[:]
   groups.extend(sd.groups)
   for group in groups:
      if group.name not in names:
         names.add(group.name)
         continue

      i = 1
      group.unique_name = group.name
      while group.unique_name in names:
         group.unique_name = group.name + " %d" % i
         i += 1
      names.add(group.unique_name)

   # header file
   WriteC4DtoAHeaderFile(sd, name, build_dir)
   # resource file
   WriteC4DtoAResFile(sd, name, build_dir)
   # string file
   WriteC4DtoAStringFile(sd, name, build_dir)


#make list of groups and tabs
def buildGroupList(sd, el):   
   if isinstance(el, Group):
      if not isinstance(el, Tab):      
         sd.groups.append(el)

   if isinstance(el, Tab):    
      sd.tabs.append(el)

   if el.children:
      for e in el.children:
         buildGroupList(sd, e)
   

#find the total number of child groups and parameters of a Tab
def findTotalChildGroups(tab):
   if tab.children:
      totalchildren=len(tab.children)  
   
      for t in tab.children:
         totalchildren+=findTotalChildGroups(t)

      return totalchildren
   else:
      return 0 


#print ordered list of all parameters and groups   
def writeParameterOrder(f, grp, first=1, orderc=1):
   firstfolder=first
   ordercount=orderc
   for gr in grp.children:
      if isinstance(gr, Parameter):
         f.write('%s ' % gr.name)

      if isinstance(gr, Group) and not isinstance(gr, Tab):
         f.write('%s ' % gr.ident)
         writeParameterOrder(f,gr,firstfolder,ordercount)
      
      if isinstance(gr, Tab):
         if firstfolder == 1:
            f.write('folder1 ')
            firstfolder=0
         ordercount+=1
         f.write('"\n\thoudini.order%d STRING "' % ordercount)
      
         writeParameterOrder(f,gr, firstfolder,ordercount)  

class HoudiniEntry:
   name = ''
   etype = ''
   def __init__(self, nm, et):
      self.name = nm
      self.etype = et

class HoudiniFolder:
   name = ''
   group_list = []

def WalkHoudiniHeader(grp, group_list):
   for el in grp.children:
      if isinstance(el, Group):
         group_list.append(HoudiniEntry(el.name, 'group'))
         WalkHoudiniHeader(el, group_list)
      elif isinstance(el, Parameter):
         group_list.append(HoudiniEntry(el.name, 'param'))

def WriteHoudiniHeader(sd, f):
   root_order_list = []
   folder_ROOT_list = []
   for el in sd.root.children:
      if isinstance(el, Group):
         new_folder = HoudiniFolder()
         new_folder.name = el.name
         group_list = []
         WalkHoudiniHeader(el, group_list)
         new_folder.group_list = group_list
         folder_ROOT_list.append(new_folder)
      elif isinstance(el, Parameter):
         root_order_list.append(HoudiniEntry(el.name, 'param'))

   # tell houdini how many groups we have and how big they are. in advance. because houdini can be dumb too sometimes
   folder_ROOT_list_STR = 'houdini.parm.folder.ROOT STRING "'
   for folder in folder_ROOT_list:
      folder_ROOT_list_STR += '%s;%d;' % (folder.name, len(folder.group_list))
   folder_ROOT_list_STR += '"'
   writei(f, folder_ROOT_list_STR, 1)

   # now tell houdini about all the headings we have...
   heading_list = []
   heading_idx = 0
   for folder in folder_ROOT_list:
      for el in folder.group_list:
         if el.etype == 'group':
            header_str = 'houdini.parm.heading.h%d STRING "%s"' % (heading_idx, el.name)
            heading_idx += 1
            writei(f, header_str, 1)

   # write main houdini ordering
   order = 'houdini.order STRING "'
   for entry in root_order_list:
      order += entry.name
      order += ' '
   order += ' ROOT"'

   writei(f, order, 1)

   # now tell houdini what's actually in the damn groups
   order_idx = 2
   heading_idx = 0
   for folder in folder_ROOT_list:
      order_str = 'houdini.order%d STRING "' % order_idx
      for el in folder.group_list:
         if el.etype == 'group':
            order_str = '%s h%d' % (order_str, heading_idx)
            heading_idx += 1
         else:
            order_str = '%s %s' % (order_str, el.name)
      order_str += '"'
      writei(f, order_str, 1)
      order_idx += 1



def WriteMDTHeader(sd, f): 
   writei(f, '[node %s]' % sd.name, 0)
   writei(f, 'desc STRING "%s"' % sd.description, 1)
   if sd.c4d_classification: writei(f, 'c4d.classification STRING "%s"' % sd.c4d_classification, 1)
   if sd.c4d_menu: writei(f, 'c4d.menu STRING "%s"' % sd.c4d_menu, 1)
   if sd.c4d_command_id: writei(f, 'c4d.command_id INT %s' % sd.c4d_command_id, 1)
   writei(f, 'maya.name STRING "%s"' % sd.name, 1)
   writei(f, 'maya.classification STRING "%s"' % sd.maya_classification, 1)
   writei(f, 'maya.id INT %s' % sd.maya_id, 1)
   #writei(f, 'houdini.label STRING "%s"' % sd.name, 1)
   if sd.houdini_icon is not None:
      writei(f, 'houdini.icon STRING "%s"' % sd.houdini_icon, 1)
   writei(f, 'houdini.category STRING "alShaders"' ,1)
   writei(f, 'houdini.help_url STRING "%s"' % sd.help_url ,1)

   WriteHoudiniHeader(sd, f)

   
def WriteMTDParam(f, name, ptype, value, d):
   if value == None:
      return

   if ptype == 'bool':
      if value:
         bval = "TRUE"
      else:
         bval = "FALSE"
      writei(f, '%s BOOL %s' % (name, bval), d)
   elif ptype == 'int':
      writei(f, '%s INT %d' % (name, value), d)
   elif ptype == 'float':
      writei(f, '%s FLOAT %r' % (name, value), d)
   elif ptype == 'string':
      writei(f, '%s STRING "%s"' % (name, value), d)

def WriteMTD(sd, fn):      
   
   #build tab and group lists
   for e in sd.root.children:
      buildGroupList(sd, e)
   
   #make all groupident unique
   sd.uniqueGroupIdents()  

   f = open(fn, 'w') 
   
   WriteMDTHeader(sd, f)   

   writei(f, '')

   for p in sd.parameters:
      writei(f, '[attr %s]' % p.name, 1)
      writei(f, 'houdini.label STRING "%s"' % p.label, 2)
      WriteMTDParam(f, "min", "float", p.mn, 2)
      WriteMTDParam(f, "max", "float", p.mx, 2)
      WriteMTDParam(f, "softmin", "float", p.smn, 2)
      WriteMTDParam(f, "softmax", "float", p.smx, 2)               
      WriteMTDParam(f, "desc", "string", p.description, 2)
      WriteMTDParam(f, "linkable", "bool", p.connectible, 2)

   for a in sd.aovs:
      writei(f, '[attr %s]' % a.name, 1)
      writei(f, 'houdini.label STRING "%s"' % a.label, 2)
      if (a.ptype == 'rgb'):
         writei(f, 'aov.type INT 0x05', 2)
      elif (a.ptype == 'rgba'):
         writei(f, 'aov.type INT 0x06', 2)
      writei(f, 'aov.enable_composition BOOL TRUE', 2)               
      writei(f, 'default STRING "%s"' % a.name[4:], 2)

   f.close()

def WalkArgs(el, f, d):
   if isinstance(el, Group):
      writei(f, '<page name="%s">' % (el.name), d)

      if el.children:
         for e in el.children:
            WalkArgs(e, f, d+1)

      writei(f, '</page>', d)
   elif isinstance(el, AOV):
      writei(f, '<param name="%s" label="%s"/>' % (el.name, el.label), d)
   elif isinstance(el, Parameter):
      if el.ptype == 'bool':
         writei(f, '<param name="%s" label="%s" widget="checkBox"/>' % (el.name, el.label), d)
      elif el.ptype == 'int':
         writei(f, '<param name="%s" label="%s" int="True"/>' % (el.name, el.label), d)
      elif el.ptype == 'rgb':
         writei(f, '<param name="%s" label="%s" widget="color"/>' % (el.name, el.label), d)
      elif el.ptype == 'enum':
         writei(f, '<param name="%s" label="%s" widget="popup">' % (el.name, el.label), d)
         
         writei(f, '<hintlist name="options">', d+1)
         
         for en in el.enum_names:
            writei(f, '<string value="%s"/>' % en, d+2)  

         writei(f, '</hintlist>', d+1)
         writei(f, '</param>', d)
      else:
         writei(f, '<param name="%s" label="%s"/>' % (el.name, el.label), d)


def WriteArgs(sd, fn):
   f = open(fn, 'w')
   writei(f, '<args format="1.0">', 0)

   for e in sd.root.children:
      WalkArgs(e, f, 0)

   writei(f, '</args>', 0)

def getSPDLTypeName(t):
   if t == 'bool':
      return 'boolean'
   elif t == 'int':
      return 'integer'
   elif t == 'float':
      return 'scalar'
   elif t == 'rgb':
      return 'color'
   elif t == 'enum':
      return 'string'
   else:
      return t

def WriteSPDLParameter(f, p):
   writei(f, 'Parameter "%s" input' % p.name, 1)
   writei(f, '{', 1)
   writei(f, 'GUID = "{%s}";' % uuid.uuid4(), 2)
   if isinstance(p, AOV):
      writei(f, 'Type = string;', 2)
      si_aov_name = p.name.replace('aov_', 'arnold_').title()
      writei(f, 'Value = "%s";' % si_aov_name, 2)
   else:
      writei(f, 'Type = %s;' % getSPDLTypeName(p.ptype), 2)
      if p.connectible:
         writei(f, 'Texturable = on;', 2)
      else:
         writei(f, 'Texturable = off;', 2)
      if p.ptype == 'bool':
         writei(f, 'Value = %s;' % str(p.default).lower(), 2)
      elif p.ptype == 'int':
         writei(f, 'Value = %d;' % p.default, 2)
      elif p.ptype == 'float':
         writei(f, 'Value = %f;' % p.default, 2)
      elif p.ptype == 'rgb' or p.ptype == 'vector':
         writei(f, 'Value = %f %f %f;' % p.default, 2)
      elif p.ptype == 'string':
         writei(f, 'Value = "%s";' % p.default, 2)

   writei(f, '}', 1)

def WriteSPDLPropertySet(f, sd):
   writei(f, 'PropertySet "%s_pset"' % sd.name)
   writei(f, '{')


   writei(f, 'Parameter "out" output', 1)
   writei(f, '{', 1)
   writei(f, 'GUID = "{%s}";' % uuid.uuid4(), 2)
   writei(f, 'Type = %s;' % getSPDLTypeName(sd.output), 2)
   writei(f, '}', 1)

   for p in sd.parameters:
      WriteSPDLParameter(f, p)

   for a in sd.aovs:
      WriteSPDLParameter(f, a)

   writei(f, '}')

def writeSPDLDefault(f, p, i):
   writei(f, '%s' % p.name, 1)
   writei(f, '{', 1)

   writei(f, 'Name = "%s";' % p.label, 2)
   if p.ptype == 'rgb':
      writei(f, 'UIType = "rgb";', 2)
   elif p.ptype=='enum':
      writei(f, 'UIType = "Combo";', 2)
      writei(f, 'Items', 2)
      writei(f, '{', 2)
      for ev in p.enum_names:
         writei(f, '"%s" = "%s";' % (ev, ev), 3)
      writei(f, '}', 2)
   elif p.smn and p.smx:
      writei(f, 'UIRange = %f To %f;' % (p.smn, p.smx), 2)


   if p.connectible:
      writei(f, 'Commands = "{F5C75F11-2F05-11D3-AA95-00AA0068D2C0}";', 2)

   writei(f, '}', 1)

def WalkSPDLDefault(f, el, d):

   if isinstance(el, Group):
      if isinstance(el, Tab):
         writei(f, 'Tab "%s"' % el.name, d)
      else:
         writei(f, 'Group "%s"' % el.name, d)
      writei(f, '{', d)
      if el.children:
         for e in el.children:
            WalkSPDLDefault(f, e, d+1)
      writei(f, '}', d)

   elif isinstance(el, Parameter):
      writei(f, '%s;' % el.name, d)

def WalkSPDLTestChildren(el):
   childrenHaveMembers = False
   if el.children:
      for e in el.children:
         if isinstance(e, Group):
            childrenHaveMembers = childrenHaveMembers or WalkSPDLTestChildren(e)
         elif isinstance(e, Parameter):
            childrenHaveMembers = childrenHaveMembers or e.connectible
   return childrenHaveMembers

def WalkSPDLRender(f, el, d):
   if isinstance(el, Group) and WalkSPDLTestChildren(el):
      writei(f, 'Group "%s"' % el.name, d)
      writei(f, '{', d)
      if el.children:
         for e in el.children:
            WalkSPDLRender(f, e, d+1)
      writei(f, '}', d)

   elif isinstance(el, Parameter) and el.connectible:
      writei(f, '%s;' % el.name, d)
      

def WriteSPDL(sd, fn):
   f = open(fn, 'w')
   writei(f, 'SPDL')
   writei(f, 'Version = "2.0.0.0";')
   writei(f, 'Reference = "{%s}";' % uuid.uuid4())

   WriteSPDLPropertySet(f, sd)

   writei(f, 'MetaShader "%s_meta"' % sd.name)
   writei(f, '{')
   writei(f, 'Name = "%s";' % sd.soft_name, 1)
   writei(f, 'Type = %s;' % sd.soft_classification, 1)
   writei(f, 'Renderer "mental ray"', 1)
   writei(f, '{', 1)
   writei(f, 'Name = "%s";' % sd.soft_name, 2)
   writei(f, 'Filename = "%s";' % sd.name, 2)
   writei(f, 'Options', 2)
   writei(f, '{', 2)
   writei(f, '"version" = 1;', 3)
   writei(f, '}', 2)
   writei(f, '}', 1)
   writei(f, '}')

   # Begin Defaults
   writei(f, 'Defaults')
   writei(f, '{')

   cmd = uuid.uuid4()

   for p in sd.parameters:
      writeSPDLDefault(f, p, cmd)

   for a in sd.aovs:
      writei(f, '%s' % a.name, 1)
      writei(f, '{', 1)
      writei(f, 'Name = "%s";' % a.name, 2)
      writei(f, '}', 1)

   writei(f, '}')
   # End Defaults

   #Begin Layout
   writei(f, 'Layout "Default"')
   writei(f, '{')
   
   for el in sd.root.children:
      WalkSPDLDefault(f, el, 1)

   writei(f, '}')
   writei(f, 'Layout "RenderTree"')
   writei(f, '{')
   
   for el in sd.root.children:
      WalkSPDLRender(f, el, 1)

   writei(f, '}')
   # End Layout

   writei(f, 'Plugin = Shader')
   writei(f, '{')
   writei(f, 'Filename = "%s";' % sd.name, 1)
   writei(f, '}')

   f.close()

def WriteHTMLParam(f, el, d):
   el_link="link='false'"
   if el.connectible:
      el_link=''
   el_range=''
   if el.mn != None and el.mx != None:
      el_range="range='[%s, %s]'" % (str(el.mn), str(el.mx))
   el_default="default='%s'" % str(el.default)
   if el.ptype == 'rgb' and not isinstance(el, AOV):
      el_default="default='rgb(%d, %d, %d)'" % (int(pow(el.default[0],1/2.2) * 255), int(pow(el.default[1],1/2.2) * 255), int(pow(el.default[2],1/2.2) * 255))
   el_fig=""
   if el.fig != None and el.figc != None:
      el_fig = "fig='%s' figc='%s'" % (el.fig, el.figc)

   el_presets=""
   if el.presets is not None:
      if el.ptype == 'float':
         el_presets = "presets='"
         first = True
         for k in sorted(el.presets, key=el.presets.get):
            if first:
               first = False
            else:
               el_presets += "|"
               
            el_presets += "%s:%.3f" % (k, el.presets[k])
         el_presets += "'"
      elif el.ptype == 'rgb':
         el_presets = "presets='"
         first = True
         for k in sorted(el.presets):
            if first:
               first = False
            else:
               el_presets += "|"
            el_presets += "%s:(%.5f, %.5f, %.5f)" % (k, el.presets[k][0], el.presets[k][1], el.presets[k][2])
         el_presets += "'"

   writei(f, "<div class='param' name='%s' label='%s' type='%s'%s desc='%s' %s %s %s %s></div>"
            % (el.name, el.label, el.ptype, el_default, el.description, el_range, el_link, el_fig, el_presets), d)

def WalkHTML(f, el, d):
   if isinstance(el, Group):
      writei(f, "<h3 class='section'><span>%s</span></h3>" % el.name, d)
      writei(f, "<div class='section-desc'>%s</div>" % el.description, d)
      if el.children:
         for e in el.children:
            WalkHTML(f, e, d+1)

   elif isinstance(el, Parameter):
      WriteHTMLParam(f, el, d)

def WalkHTMLRoot(f, el, d):
   if isinstance(el, Group):
      expand_class='expander-trigger'
      if el.collapse:
         expand_class='expander-trigger expander-hidden'
      writei(f, "<div class='expander'>", d)
      writei(f, "<div class='%s'>%s</div>" % (expand_class, el.name), d)
      writei(f, "<div class='expander-content'>", d)
      writei(f, "<div class='section-desc'>%s</div>" % el.description, d)

      if el.children:
         for e in el.children:
            WalkHTML(f, e, d+1)

      writei(f, "</div> <!-- end expander-content %s -->" % el.name, d)
      writei(f, "</div> <!-- end expander %s -->" % el.name, d)

   elif isinstance(el, Parameter):
      WriteHTMLParam(f, el, d)

def WriteHTML(sd, fn):
   f = open(fn, 'w')

   f.write(
      """
<!DOCTYPE html>
<html>
<head>
   <meta charset='utf-8'>
   <title>%s Reference Guide</title>
   <link rel="stylesheet" href="../css/app.css">
   <script src="http://ajax.googleapis.com/ajax/libs/jquery/1.11.2/jquery.min.js"></script>
   <script src="../js/navbar.js"></script>
   <script src="../js/tabs.js"></script>
   <script src="../js/param.js"></script>
</head>

<header class="navigationnofix fixed" role="banner">
  <div class="navigation-wrapper">
    <a href="index.html" class="logo">
      <img src="img/als-logo-white.svg" alt="Logo Image">
    </a>
    <a href="javascript:void(0)" class="navigation-menu-button" id="js-mobile-menu"><i class="fa fa-bars fa-2x"></i></a>
    <nav role="navigation">
      <ul id="js-navigation-menu" class="navigation-menu show">
        <li class="nav-link"><a href="index.html">Home</a></li>
            <li class="nav-link"><a href="javascript:void(0)">Tutorials</a></li>
            <li class="nav-link more"><a href="javascript:void(0)">Surfaces</a>
            <ul class="submenu">
               <li><a href="alSurface.html">alSurface</a></li>
               <li><a href="alLayer.html">alLayer</a></li>
               <li><a href="alHair.html">alHair</a></li>
               <li><a href="alCel.html">alCel</a></li>
            </ul>
            </li>
            <li class="nav-link more"><a href="javascript:void(0)">Textures</a>
            <ul class="submenu">
               <li><a href="alFractal.html">alFractal</a></li>
               <li><a href="alCellNoise.html">alCellNoise</a></li>
               <li><a href="alGaborNoise.html">alGaborNoise</a></li>
               <li><a href="alFlowNoise.html">alFlowNoise</a></li>
               <li><a href="alPattern.html">alPattern</a></li>
               <li><a href="alTriplanar.html">alTriplanar</a></li>
               <li><a href="alBlackbody.html">alBlackbody</a></li>
            </ul>
            </li>
            <li class="nav-link more"><a href="javascript:void(0)">Utilities</a>
            <ul class="submenu">
               <li><a href="alCombineFloat.html">alCombineFloat</a></li>
               <li><a href="alCombineColor.html">alCombineColor</a></li>
               <li><a href="alRemapFloat.html">alRemapFloat</a></li>
               <li><a href="alRemapColor.html">alRemapColor</a></li>
               <li><a href="alLayerFloat.html">alLayerFloat</a></li>
               <li><a href="alLayerColor.html">alLayerColor</a></li>
               <li><a href="alSwitchFloat.html">alSwitchFloat</a></li>
               <li><a href="alSwitchColor.html">alSwitchColor</a></li>
               <li><a href="alCache.html">alCache</a></li>
               <li><a href="alCurvature.html">alCurvature</a></li>
               <li><a href="alInputVector.html">alInputVector</a></li>
               <li><a href="alInputScalar.html">alInputScalar</a></li>
               <li><a href="alJitterColor.html">alJitterColor</a></li>
      </ul>
    </nav>
    
  </div>
</header>

<body>
   <div class='shader-content'>
      <div class='article-block'>
         <article class='type-system-sans-dark'>
            <p class='type'>%s</p>
            <h2>%s</h2>
            <p>%s This document is a reference guide to the parameters.</p>
         </article>
      </div> <!-- end article-block -->
   
      <div class='parameter-block'>
      <div class='parameter-content'>
"""
   % (sd.name, sd.name, sd.intro, sd.description)
   )

   # parameters
   for el in sd.root.children:
      WalkHTMLRoot(f, el, 2)

   f.write(
"""
      </div> <!-- end parameter-content -->
      </div> <!-- end parameter-block -->
   </div> <!-- end shader-content -->
</body>
"""
)

   f.close()

def remapControls(sd):
   with group(sd, 'Remap', collapse=True, description='These controls allow you to remap the shader result.'):
      sd.parameter('RMPinputMin', 'float', 0.0, label='Input min', description='Sets the minimum input value. Use this to pull values outside of 0-1 into a 0-1 range.')
      sd.parameter('RMPinputMax', 'float', 1.0, label='Input max', description='Sets the maximum input value. Use this to pull values outside of 0-1 into a 0-1 range.')

      with group(sd, 'Contrast', collapse=False):
         sd.parameter('RMPcontrast', 'float', 1.0, label='Contrast', description='Scales the contrast of the input signal.')
         sd.parameter('RMPcontrastPivot', 'float', 0.18, label='Pivot', description='Sets the pivot point around which the input signal is contrasted.')

      with group(sd, 'Bias and gain', collapse=False):
         sd.parameter('RMPbias', 'float', 0.5, label='Bias', description='Bias the signal higher or lower. Values less than 0.5 push the average lower, values higher than 0.5 push it higher.')
         sd.parameter('RMPgain', 'float', 0.5, label='Gain', description='Adds gain to the signal, in effect a different form of contrast. Values less than 0.5 increase the gain, values greater than 0.5 decrease it.')

      sd.parameter('RMPoutputMin', 'float', 0.0, label='Output min', description='Sets the minimum value of the output. Use this to scale a 0-1 signal to a new range.')
      sd.parameter('RMPoutputMax', 'float', 1.0, label='Output max', description='Sets the maximum value of the output. Use this to scale a 0-1 signal to a new range.')

      with group(sd, 'Clamp', collapse=False):
         sd.parameter('RMPclampEnable', 'bool', False, label='Enable', description='When enabled, will clamp the output to Min-Max.')
         sd.parameter('RMPthreshold', 'bool', False, label='Expand', description='When enabled, will expand the clamped range to 0-1 after clamping.')
         sd.parameter('RMPclampMin', 'float', 0.0, label='Min', description='Minimum value to clamp to.')
         sd.parameter('RMPclampMax', 'float', 1.0, label='Max', description='Maximum value to clamp to.')

# Main. Load the UI file and build UI templates from the returned structure
if __name__ == '__main__':
   if len(sys.argv) < 7:
      print 'ERROR: must supply exactly ui source input and mtd, ae, spdl and args outputs'
      sys.exit(1)

   ui = ShaderDef()
   globals_dict = {'ui':ui}
   execfile(sys.argv[1], globals_dict)

   if not isinstance(ui, ShaderDef):
      print 'ERROR: ui object is not a ShaderDef. Did you assign something else to it by mistake?'
      sys.exit(2)

   WriteMTD(ui, sys.argv[2])  
   WriteAETemplate(ui, sys.argv[3])
   WriteSPDL(ui, sys.argv[4])
   WriteArgs(ui, sys.argv[5])

   # C4DtoA resource files
   name = os.path.basename(os.path.splitext(sys.argv[1])[0])
   build_dir = sys.argv[6] if len(sys.argv) > 6 else os.path.abspath("")
   WriteC4DtoAResourceFiles(ui, name, build_dir)

   # HTML
   WriteHTML(ui, sys.argv[7])


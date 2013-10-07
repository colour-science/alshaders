#! /usr/bin/env python

import sys
import uuid

def enum(*sequential, **named):
	enums = dict(zip(sequential, range(len(sequential))), **named)
	reverse = dict((value, key) for key, value in enums.iteritems())
	enums['reverse_mapping'] = reverse
	return type('Enum', (), enums)

class UiElement:
	name = ''
	description = None
	parent = None

class Group(UiElement):
	collapse = True
	children = None

	def __init__(self, name, collapse=True, description=None):
		self.name = name
		self.collapse = collapse
		self.description = description

	def __str__(self):
		return 'Group %s' % self.name


class Tab(Group):
	def __init__(self, name, collapse=True, description=None):
		self.name = name
		self.collapse = collapse
		self.description = description

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

	def __init__(self, name, ptype, label=None, description=None, mn=None, mx=None, smn=None, smx=None, connectible=True, enum_names=None):
		self.name = name
		self.ptype = ptype
		self.description = description
		if not label:
			label = name
		self.label = label
		self.mn = mn
		self.mx = mx
		self.smn = smn
		self.smx = smx
		if ptype == 'enum':
			self.enum_names = enum_names

		if ptype == 'bool' or ptype == 'string' or ptype == 'int' or ptype == 'enum':
			self.connectible = False
		else:
			self.connectible = connectible

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
	description = None
	output = None
	maya_name = None
	maya_classification = None
	maya_id = None
	maya_swatch = False
	maya_bump = False
	maya_matte = False
	soft_name = None
	soft_classification = None
	soft_version = 1
	hierarchy_depth = 0
	root = None
	current_parent = None
	parameters = []
	aovs = []

	def __init__(self):
		self.root = Group('__ROOT__', False)
		self.current_parent = self.root

	def setMayaClassification(self, mc):
		self.maya_classification = mc

	def setDescription(self, d):
		self.description = d

	def setOutput(self, o):
		self.output = o

	def beginGroup(self, name, collapse=True, description=''):
		g = None
		if self.current_parent is self.root:
			g = Tab(name, collapse, description)
		else:
			g = Group(name, collapse, description)

		if not self.current_parent.children:
			self.current_parent.children = [g]
		else:
			self.current_parent.children.append(g)

		g.parent = self.current_parent
		self.current_parent = g

	def endGroup(self):
		self.current_parent = self.current_parent.parent

	def parameter(self, name, ptype, label=None, description=None, mn=None, mx=None, smn=None, smx=None, connectible=None, enum_names=None):
		p = Parameter(name, ptype, label, description, mn, mx, smn, smx, connectible, enum_names)
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
		self.description = d['description']
		self.output = d['output']
		self.maya_name = d['maya_name']
		self.maya_classification = d['maya_classification']
		self.maya_id = d['maya_id']
		self.maya_swatch = d['maya_swatch']
		self.maya_bump = d['maya_bump']
		self.maya_matte = d['maya_matte']
		self.soft_name = d['soft_name']
		self.soft_classification = d['soft_classification']
		self.soft_version = d['soft_version']

class group:
	name = ''
	collapse = True
	description = ''
	ui = None

	def __init__(self, ui, name, collapse=True, description=None):
		self.ui = ui
		self.name = name
		self.collapse = collapse
		self.description = description

	def __enter__(self):
		self.ui.beginGroup(self.name, self.collapse, self.description)

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
		writei(f, 'self.addControl("%s", label="%s")' % (el.name, el.label), d)

def WriteAETemplate(sd, fn):
	f = open(fn, 'w')
	writei(f, 'import pymel.core as pm', 0)
	writei(f, 'from alShaders import alShadersTemplate\n', 0)

	writei(f, 'class AE%sTemplate(alShadersTemplate):' % sd.name, 0)
	writei(f, 'def setup(self):', 1)

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
	writei(f, 'self.endScrollLayout()', 2)	#end main scrollLayout

	f.close()

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
		writei(f, '%s STRING "%r"' % (name, value), d)

def WriteMTD(sd, fn):
	f = open(fn, 'w')
	writei(f, '[node %s]' % sd.name, 0)
	writei(f, 'desc STRING "%s"' % sd.description, 1)
	writei(f, 'maya.name STRING "%s"' % sd.name, 1)
	writei(f, 'maya.classification STRING "%s"' % sd.maya_classification, 1)
	writei(f, 'maya.id INT %s' % sd.maya_id, 1)
	writei(f, '')

	for p in sd.parameters:
		writei(f, '[attr %s]' % p.name, 1)
		WriteMTDParam(f, "min", "float", p.mn, 2)
		WriteMTDParam(f, "max", "float", p.mx, 2)
		WriteMTDParam(f, "softmin", "float", p.smn, 2)
		WriteMTDParam(f, "softmax", "float", p.smx, 2)
		WriteMTDParam(f, "desc", "string", p.description, 2)
		WriteMTDParam(f, "linkable", "bool", p.connectible, 2)

	for a in sd.aovs:
		writei(f, '[attr %s]' % a.name, 1)
		writei(f, 'aov.type INT 0x05', 2)
		writei(f, 'aov.enable_composition BOOL TRUE', 2)
		writei(f, 'default STRING "%s"' % a.name[4:], 2)

	f.close()

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
		writei(f, 'Commands = "{F5C75F11-2F05-11D3-AA95-00AA0068D2C0}";' % i, 2)

	writei(f, '}', 1)

def WalkSPDL(f, el, d, wc):

	if isinstance(el, Group):
		if isinstance(el, Tab) and wc:
			writei(f, 'Tab "%s"' % el.name, d)
		else:
			writei(f, 'Group "%s"' % el.name, d)
		writei(f, '{', d)
		if el.children:
			for e in el.children:
				WalkSPDL(f, e, d+1, wc)
		writei(f, '}', d)

	elif isinstance(el, Parameter):
		if wc or el.connectible:
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
		WalkSPDL(f, el, 1, True)

	writei(f, '}')
	writei(f, 'Layout "RenderTree"')
	writei(f, '{')
	
	for el in sd.root.children:
		WalkSPDL(f, el, 1, False)

	writei(f, '}')
	# End Layout

# 	#Begin Logic
# 	writei(f, 'Logic')
# 	writei(f, '{')
	
# 	if sd.aovs:
# 		f.write("""
# 	Sub OnInit
# 		Dim oRenderChannels
# 		Set oRenderChannels = ActiveProject.ActiveScene.PassContainer.Properties( "Scene Render Options" ).RenderChannels

# 		If TypeName(oRenderChannels) = "Nothing" Then
# 			LogMessage "Scene Render Options property not found. Can't enumerate render channels.", siError
# 			Exit Sub
# 		End If
# 		If oRenderChannels.Count = 0 Then
# 			LogMessage "No render channels defined.", siError
# 			Exit Sub
# 		End If

# 		Dim		idx
# 		ReDim oChannelList( oRenderChannels.Count * 2 + 1 )

# 		idx = 0
# 		for each oChannel in oRenderChannels
# 			If oChannel.ChannelType = siRenderChannelColorType And oChannel.UserDefined Then
# 				oChannelList( idx * 2 + 0 ) = oChannel.Name
# 				oChannelList( idx * 2 + 1 ) = oChannel.Name
# 				idx = idx + 1
# 			End If
# 		next

# 		ReDim Preserve oChannelList( idx * 2 - 1 )

# 		Dim oChannelCombo
# """)
# 		for a in sd.aovs:
# 			writei(f, 'Set oChannelCombo = PPG.PPGLayout.Item("%s")' % a.name, 2)
# 			writei(f, 'oChannelCombo.UIItems = oChannelList', 2)

# 		f.write("""

# 		PPG.Refresh
# 	End Sub
# """)
# 	writei(f, '}')
# 	#End Logic

	writei(f, 'Plugin = Shader')
	writei(f, '{')
	writei(f, 'Filename = "%s";' % sd.name, 1)
	writei(f, '}')

	f.close()

def remapControls(sd):
	with group(sd, 'Remap', collapse=True):
		sd.parameter('RMPinputMin', 'float', label='Input min')
		sd.parameter('RMPinputMax', 'float', label='Input max')

		with group(sd, 'Contrast', collapse=False):
			sd.parameter('RMPcontrast', 'float', label='Contrast')
			sd.parameter('RMPcontrastPivot', 'float', label='Pivot')

		with group(sd, 'Bias and gain', collapse=False):
			sd.parameter('RMPbias', 'float', label='Bias')
			sd.parameter('RMPgain', 'float', label='Gain')

		sd.parameter('RMPoutputMin', 'float', label='Output min')
		sd.parameter('RMPoutputMax', 'float', label='Output max')

		with group(sd, 'Clamp', collapse=False):
			sd.parameter('RMPclampEnable', 'float', label='Enable')
			sd.parameter('RMPthreshold', 'float', label='Expand')
			sd.parameter('RMPclampMin', 'float', label='Min')
			sd.parameter('RMPclampMax', 'float', label='Max')

# Main. Load the UI file and build UI templates from the returned structure
if __name__ == '__main__':
	if len(sys.argv) != 5:
		print 'ERROR: must supply exactly ui source input and mtd, ae and spdl outputs'
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

	



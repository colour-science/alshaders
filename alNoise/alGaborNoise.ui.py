import uigen

ui.shader({
	'name':'alGaborNoise',
	'description':'Gabor noise pattern generator',
	'output':'rgb',
	'maya_name':'alGaborNoise',
	'maya_classification':'texture',
	'maya_id':'0x00116410',
	'maya_swatch':True,
	'maya_matte':False,
	'maya_bump':False,
	'soft_name':'ALS_GaborNoise',
	'soft_classification':'texture',
	'soft_version':1
})

ui.parameter('space', 'enum', 'Space', enum_names=[
	"world",
	"object",
	"Pref",
	"UV"
])
ui.parameter('frequency', 'float', 'Frequency', connectible=False)
ui.parameter('octaves', 'int', 'Octaves', connectible=False)
ui.parameter('anisotropy', 'enum', 'Anisotropy', enum_names=[
	"isotropic",
	"anisotropic",
	"hybrid"
])
ui.parameter('anisotropyDirection', 'vector', 'Anisotropy direction')
ui.parameter('bandwidth', 'float', 'Bandwidth', connectible=False)
ui.parameter('angle', 'float', 'Angle', connectible=False)
ui.parameter('advection', 'float', 'Advection', connectible=False)
ui.parameter('turbulent', 'bool', 'Turbulent')

uigen.remapControls(ui)

ui.parameter('color1', 'rgb', 'Color 1')
ui.parameter('color2', 'rgb', 'Color 2')
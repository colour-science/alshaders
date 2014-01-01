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

ui.parameter('space', 'enum','world' 'Space', enum_names=[
	"world",
	"object",
	"Pref",
	"UV"
])
ui.parameter('frequency', 'float', 1.0, 'Frequency', connectible=False)
ui.parameter('anisotropy', 'enum','isotropic', 'Anisotropy', enum_names=[
	"isotropic",
	"anisotropic",
	"hybrid"
])
ui.parameter('anisotropyDirection', 'vector', (0.0, 1.0, 0.0), 'Anisotropy direction')
ui.parameter('filter', 'bool', False, 'Filter')
ui.parameter('bandwidth', 'float', 1.0, 'Bandwidth', connectible=False)
ui.parameter('impulses', 'float', 8.0, 'Bandwidth', connectible=False)
ui.parameter('turbulent', 'bool', False, 'Turbulent')

uigen.remapControls(ui)

ui.parameter('color1', 'rgb', (0.0, 0.0, 0.0), 'Color 1')
ui.parameter('color2', 'rgb', (0.0, 0.0, 0.0), 'Color 2')

ui.parameter('P', 'vector', (0.0, 0.0, 0.0), 'P')
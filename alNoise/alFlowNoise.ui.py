import uigen

ui.shader({
	'name':'alFlowNoise',
	'description':'Flow noise pattern generator',
	'output':'rgb',
	'maya_name':'alFlowNoise',
	'maya_classification':'texture',
	'maya_id':'0x0011640F',
	'maya_swatch':True,
	'maya_matte':False,
	'maya_bump':False,
	'soft_name':'ALS_FlowNoise',
	'soft_classification':'texture',
	'soft_version':1
})

ui.parameter('space', 'enum', 'world', 'Space', enum_names=[
	"world",
	"object",
	"Pref",
	"UV"
])
ui.parameter('frequency', 'float', 1.0, 'Frequency', connectible=False)
ui.parameter('octaves', 'int', 4, 'Octaves', connectible=False)
ui.parameter('lacunarity', 'float', 2.172, 'Lacunarity', connectible=False)
ui.parameter('gain', 'float', 0.5, 'Gain', connectible=False)
ui.parameter('angle', 'float', 0.0, 'Angle', connectible=False)
ui.parameter('advection', 'float', 0.25, 'Advection', connectible=False)
ui.parameter('turbulent', 'bool', False, 'Turbulent')

uigen.remapControls(ui)

ui.parameter('color1', 'rgb', (0.0, 0.0, 0.0), 'Color 1')
ui.parameter('color2', 'rgb', (0.0, 0.0, 0.0), 'Color 2')

ui.parameter('P', 'vector', (0.0, 0.0, 0.0), 'P')
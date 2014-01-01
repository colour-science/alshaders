import uigen

ui.shader({
	'name':'alFractal',
	'description':'Fractal noise pattern generator',
	'output':'rgb',
	'maya_name':'alFractal',
	'maya_classification':'texture',
	'maya_id':'0x00116409',
	'maya_swatch':True,
	'maya_matte':False,
	'maya_bump':False,
	'soft_name':'ALS_Fractal',
	'soft_classification':'texture',
	'soft_version':1
})

ui.parameter('mode', 'enum', 'scalar', 'Mode', enum_names=[
	"scalar",
	"vector"
])
ui.parameter('space', 'enum', 'world', 'Space', enum_names=[
	"world",
	"object",
	"Pref",
	"UV"
])
ui.parameter('scale', 'vector', (1.0, 1.0, 1.0), 'Vector', connectible=False)
ui.parameter('frequency', 'float', 1.0, 'Frequency', connectible=False)
ui.parameter('time', 'float', 0.0, 'Time', connectible=False)
ui.parameter('octaves', 'int', 8, 'Octaves', connectible=False)
ui.parameter('distortion', 'float', 0.0, 'Distortion', connectible=False)
ui.parameter('lacunarity', 'float', 2.121, 'Lacunarity', connectible=False)
ui.parameter('gain', 'float', 0.5, 'Gain', connectible=False)
ui.parameter('turbulent', 'bool', False, 'Turbulent')
ui.parameter('ridged', 'bool', False, 'Ridged')
ui.parameter('ridgeOffset', 'float', 0.0, 'Ridge offset', connectible=False)

uigen.remapControls(ui)

ui.parameter('color1', 'rgb', (0.0, 0.0, 0.0), 'Color 1')
ui.parameter('color2', 'rgb', (0.0, 0.0, 0.0), 'Color 2')

ui.parameter('P', 'vector', (0.0, 0.0, 0.0), 'P')
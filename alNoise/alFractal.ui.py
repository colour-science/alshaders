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

ui.parameter('mode', 'enum', 'Space', enum_names=[
	"scalar",
	"vector"
])
ui.parameter('space', 'enum', 'Space', enum_names=[
	"world",
	"object",
	"Pref",
	"UV"
])
ui.parameter('scale', 'vector', 'Vector', connectible=False)
ui.parameter('frequency', 'float', 'Frequency', connectible=False)
ui.parameter('time', 'float', 'Time', connectible=False)
ui.parameter('octaves', 'int', 'Octaves', connectible=False)
ui.parameter('distortion', 'float', 'Distortion', connectible=False)
ui.parameter('lacunarity', 'float', 'Lacunarity', connectible=False)
ui.parameter('gain', 'float', 'Gain', connectible=False)
ui.parameter('turbulent', 'bool', 'Turbulent')
ui.parameter('ridged', 'bool', 'Ridged')
ui.parameter('ridgeOffset', 'float', 'Ridge offset', connectible=False)

uigen.remapControls(ui)

ui.parameter('color1', 'rgb', 'Color 1')
ui.parameter('color2', 'rgb', 'Color 2')
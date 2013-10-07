import uigen

ui.shader({
	'name':'alCellNoise',
	'description':'Cell noise pattern generator',
	'output':'rgb',
	'maya_name':'alCellNoise',
	'maya_classification':'texture',
	'maya_id':'0x00116408',
	'maya_swatch':True,
	'maya_matte':False,
	'maya_bump':False,
	'soft_name':'ALS_CellNoise',
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
ui.parameter('mode', 'enum', 'Mode', enum_names=[
	"features",
	"chips"
])

with uigen.group(ui, 'Features', False):
	ui.parameter('octaves', 'int', 'Octaves', connectible=False)
	ui.parameter('lacunarity', 'int', 'Lacunarity', connectible=False)
	ui.parameter('randomness', 'int', 'Randomness', connectible=False)

	with uigen.group(ui, 'Feature weights', False):
		ui.parameter('f1w', 'float', 'F1')
		ui.parameter('f2w', 'float', 'F2')
		ui.parameter('f3w', 'float', 'F3')
		ui.parameter('f4w', 'float', 'F4')
		ui.parameter('mynkowskiShape', 'float', 'Shape')

		uigen.remapControls(ui)

		ui.parameter('color1', 'rgb', 'Color 1')
		ui.parameter('color2', 'rgb', 'Color 2')

with uigen.group(ui, 'Chips', False):
	ui.parameter('smoothChips', 'bool', 'Smooth chips')
	ui.parameter('chipColor1', 'rgb', 'Chip color 1')
	ui.parameter('chipProb1', 'float', 'Chip probability 1')
	ui.parameter('chipColor2', 'rgb', 'Chip color 2')
	ui.parameter('chipProb2', 'float', 'Chip probability 2')
	ui.parameter('chipColor3', 'rgb', 'Chip color 3')
	ui.parameter('chipProb3', 'float', 'Chip probability 3')
	ui.parameter('chipColor4', 'rgb', 'Chip color 4')
	ui.parameter('chipProb4', 'float', 'Chip probability 4')
	ui.parameter('chipColor5', 'rgb', 'Chip color 5')
	ui.parameter('chipProb5', 'float', 'Chip probability 5')

ui.parameter('P', 'vector', 'P')

import uigen

ui.shader({
	'name':'alCombineFloat',
	'description':'Combines (multiply, add, mix etc.) two floats together',
	'output':'float',
	'maya_name':'alCombineFloat',
	'maya_classification':'texture',
	'maya_id':'0x00116406',
	'maya_swatch':True,
	'maya_matte':False,
	'maya_bump':False,
	'soft_name':'ALS_CombineFloat',
	'soft_classification':'texture',
	'soft_version':1
})

ui.parameter('input1', 'float', 0.0)
ui.parameter('input2', 'float', 0.0)
ui.parameter('input3', 'float', 0.0)
ui.parameter('combineOp', 'enum','multiply 1*2', 'Combine Op', enum_names=["multiply 1*2",
	"add 1+2",
	"divide 1/2",
	"subtract 1-2",
	"lerp(1, 2, 3)"])
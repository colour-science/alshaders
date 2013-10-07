import uigen

ui.shader({
	'name':'alCombineColor',
	'description':'Combines (multiply, add, mix etc.) two colors together',
	'output':'rgb',
	'maya_name':'alCombineColor',
	'maya_classification':'texture',
	'maya_id':'0x00116405',
	'maya_swatch':True,
	'maya_matte':False,
	'maya_bump':False,
	'soft_name':'ALS_CombineColor',
	'soft_classification':'texture',
	'soft_version':1
})

ui.parameter('input1', 'rgb')
ui.parameter('input2', 'rgb')
ui.parameter('input3', 'rgb')
ui.parameter('combineOp', 'enum', 'Combine Op', enum_names=["multiply 1*2",
	"add 1+2",
	"divide 1/2",
	"subtract 1-2",
	"lerp(1, 2, 3)"])
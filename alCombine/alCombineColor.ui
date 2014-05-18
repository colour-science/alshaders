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
	'soft_version':1,
	'help_url':'https://bitbucket.org/anderslanglands/alshaders/wiki/alCombine'
})

ui.parameter('input1', 'rgb', (0.0, 0.0, 0.0))
ui.parameter('input2', 'rgb', (0.0, 0.0, 0.0))
ui.parameter('input3', 'rgb', (0.0, 0.0, 0.0))
ui.parameter('combineOp', 'enum',  'multiply 1*2', 'Combine Op', enum_names=["multiply 1*2",
	"add 1+2",
	"divide 1/2",
	"subtract 1-2",
	"lerp(1, 2, 3)"])

import uigen

ui.shader({
        'name':'alIntToFloat',
        'description':'Converts int input to float output',
        'output':'float',
        'maya_name':'alIntToFloat',
	'maya_classification':'texture',
        'maya_id':'0x0011641A',
	'maya_swatch':True,
	'maya_matte':False,
	'maya_bump':False,
        'soft_name':'ALS_IntToFloat',
	'soft_classification':'texture',
	'soft_version':1
})

ui.parameter('input', 'int', 0, 'Input')

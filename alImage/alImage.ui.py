import uigen

ui.shader({
        'name':'alImage',
        'description':'Texture node with support for hashed frame number',
	    'output':'rgb',
        'maya_name':'alImage',
	    'maya_classification':'texture',
        'maya_id':'0x00116419',
	    'maya_swatch':True,
	    'maya_matte':False,
	    'maya_bump':False,
        'soft_name':'ALS_Image',
	    'soft_classification':'texture',
	    'soft_version':1
})

ui.parameter('filename', 'string', '', 'Filename', connectible=False)
ui.parameter('frame', 'int', 0, 'Frame')

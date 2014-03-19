import uigen

ui.shader({
        'name':'alCache',
        'description':'Caches the input value for optimization of multple shader calls.',
        'output':'rgb',
        'maya_name':'alCache',
        'maya_classification':'texture',
        'maya_id':'0x00116417',
        'maya_swatch':True,
        'maya_matte':False,
        'maya_bump':False,
        'soft_name':'ALS_Cache',
        'soft_classification':'texture',
        'soft_version':1
})

ui.parameter('input', 'rgb', (0.0, 0.0, 0.0))


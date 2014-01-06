import uigen

ui.shader({
        'name':'alTriplanar',
        'description':'Projects a texture orthographically in three axis.',
        'output':'rgb',
        'maya_name':'alTriplanar',
        'maya_classification':'texture',
        'maya_id':'0x00116413',
        'maya_swatch':True,
        'maya_matte':False,
        'maya_bump':False,
        'soft_name':'ALS_Triplanar',
        'soft_classification':'texture',
        'soft_version':1
})

ui.parameter('space', 'enum', 'world', 'Space', enum_names=[
        "world",
        "object",
        "Pref"
])

ui.parameter('tiling', 'enum', 'regular', 'Tiling', enum_names=[
        "regular",
        "cellnoise"
])

ui.parameter('texture', 'string', '', 'Texture', connectible=False)

ui.parameter('blendSoftness', 'float', 0.1, 'Blend Softness', "", 0., 1., 0., 1., connectible=False)

ui.parameter('scale', 'float', 1.0, 'Scale', connectible=False)

ui.parameter('offset', 'float', 0.0, 'Offset', connectible=False)

ui.parameter('P', 'vector', (0.0, 0.0, 0.0), 'P')

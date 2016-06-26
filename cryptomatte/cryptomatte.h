#include <ai.h>
#include <string>
#include <cstring>
#include <map>
#include <ctime>
#include "MurmurHash3.h"

#define NOMINMAX // lets you keep using std::min

///////////////////////////////////////////////
//
//      Shader Data Struct
//
///////////////////////////////////////////////


struct AutomatteData {
    std::string aov_autoasset;
    std::string aov_autoobject;
    std::string aov_automaterial;
    std::string aov_cryptoasset;
    std::string aov_cryptoobject;
    std::string aov_cryptomaterial;
    AtArray * aovArray_cryptoasset;
    AtArray * aovArray_cryptoobject;
    AtArray * aovArray_cryptoMaterial;
    bool strip_mat_ns;
};

///////////////////////////////////////////////
//
//      Automatte Data
//
///////////////////////////////////////////////

#define MAX_STRING_LENGTH 255
const int MAX_CRYPTOMATTE_DEPTH = 99;

unsigned char g_pointcloud_instance_verbosity = 0;


// UINT64 djb2_hash(char *name) {
//     UINT32 len = strlen(name);
//     UINT64 hash = 5381;
//     for (UINT32 i=0; i<len; i++)
//         hash = (hash << 5) + hash + ((unsigned char) name[i]);
//     for (UINT32 i=0; i<16; i++)
//         hash = (hash << 5) + hash;
//     return hash ;
// }


bool sitoa_pointcloud_instance_handling(const char *obj_full_name, char *obj_name_out) {
    if (g_pointcloud_instance_verbosity == 0 || strstr(obj_full_name, ".SItoA.Instance.") == NULL)  {
        return false;
    }
    // AiMsgWarning("Name debug: Found a trouble name! %s", obj_full_name);
    
    char obj_name[MAX_STRING_LENGTH];
    size_t mdl_chars = strlen(obj_full_name);
    strncpy(obj_name, obj_full_name, mdl_chars);
    obj_name[ MAX_STRING_LENGTH-1 ] = '\0';

    char *instance_start = strstr(obj_name, ".SItoA.Instance.");

    char *space = strstr(instance_start, " ");
    if (space == NULL) { return false; }

    char *instance_name = &space[1];
    // AiMsgWarning("Name debug: Instance Name! %s", instance_name);
    
    char *obj_suffix2 = strstr(instance_name, ".SItoA."); 
    if (obj_suffix2 == NULL) { return false; }  
    obj_suffix2[0] = '\0';  // strip the suffix
    // AiMsgWarning("Name debug: suffix! %s", obj_suffix2);

    size_t chars_to_copy = strlen(instance_name);
    if (chars_to_copy >= MAX_STRING_LENGTH || chars_to_copy == 0) { return false; } 
    // AiMsgWarning("Name debug: Attempting the copy of %d chars! ", chars_to_copy);

    if (g_pointcloud_instance_verbosity == 2)   {
        char *frame_numbers = &instance_start[16]; // 16 chars in ".SItoA.Instance.", this gets us to the first number
        char *instance_ID = strstr(frame_numbers, ".");
        char *instance_ID_end = strstr(instance_ID, " ");
        instance_ID_end[0] = '\0';
        // AiMsgWarning("Name debug: Instance ID is %s", instance_ID);

        size_t ID_len = strlen(instance_ID);
        strncpy(&instance_name[chars_to_copy], instance_ID, ID_len);
        chars_to_copy += ID_len;
    }

    strncpy(obj_name_out, instance_name, chars_to_copy);
    obj_name_out[chars_to_copy] = '\0';
    // AiMsgWarning("Name debug: New Name! %s", obj_name_out);

    
    return true;
}


void get_clean_object_name(const char *obj_full_name, char *obj_name_out, char *mdl_name_out) { 
    char obj_name[MAX_STRING_LENGTH];
    char mdl_name[MAX_STRING_LENGTH];
    memset(obj_name, 0, MAX_STRING_LENGTH);
    memset(mdl_name, 0, MAX_STRING_LENGTH); 

    size_t mdl_chars = strlen(obj_full_name);
    if (mdl_chars > MAX_STRING_LENGTH) {
        mdl_chars = MAX_STRING_LENGTH;
    }
    
    strncpy(mdl_name, obj_full_name, mdl_chars);
    mdl_name[ MAX_STRING_LENGTH-1 ] = '\0';

    bool preempt_object_name = false;

    char *obj_postfix = strstr(mdl_name, ".SItoA.");
    if (obj_postfix != NULL) {
        // in Softimage mode
        // to do: when there are more than one way to preempt object names here, we're going to have to have some kind of loop handling that. 
        preempt_object_name = sitoa_pointcloud_instance_handling(obj_full_name, obj_name_out);
        obj_postfix[0] = '\0';
    }

    char *space_finder = strstr(mdl_name, " ");
    while (space_finder != NULL) {
        space_finder[0] = '#';
        space_finder = strstr(mdl_name, " ");
    }

    char *mdl_separator = strchr(mdl_name, ':');
    if (mdl_separator == NULL) {
        mdl_separator = strchr(mdl_name, '.');
    }

    if (mdl_separator != NULL) {
        mdl_separator[0] = '\0';
        char *obj_name_start = mdl_separator + 1;
        strncpy(obj_name, obj_name_start, strlen(obj_name_start));
    } else {
        // no namespace
        strncpy(obj_name, mdl_name, strlen(mdl_name)); // the object name is the model name in this case
        strncpy(mdl_name, "default\0", 8); // and the model name is default. 
    }

    if (!preempt_object_name) {
        strcpy(obj_name_out, obj_name);
    }
    strcpy(mdl_name_out, mdl_name);
}

void get_clean_material_name(const char *mat_full_name, char *mat_name_out, bool strip_ns) {    
    // Example: 
    //      Softimage: Sources.Materials.myLibrary_ref_library.myMaterialName.Standard_Mattes.uBasic.SITOA.25000....
    //      Maya: namespace:my_material_sg

    char mat_name[MAX_STRING_LENGTH];
    memset(mat_name, 0, MAX_STRING_LENGTH);

    if (mat_full_name != NULL) {
        size_t mat_chars = strlen(mat_full_name);
        if (mat_chars > MAX_STRING_LENGTH) {
            mat_chars = MAX_STRING_LENGTH;
        }

        strncpy(mat_name, mat_full_name, mat_chars);
        mat_name[ sizeof(mat_name)-1 ] = '\0';


        char *mat_postfix = strstr(mat_name, ".SItoA.");
        if (mat_postfix != NULL) {
            //   Sources.Materials.myLibrary_ref_library.myMaterialName.Standard_Mattes.uBasic  <<chop>> .SITOA.25000....
            mat_postfix[0] = '\0';

            char *mat_shader_name = strrchr(mat_name, '.');
            if (mat_shader_name != NULL) {
                //   Sources.Materials.myLibrary_ref_library.myMaterialName.Standard_Mattes <<chop>> .uBasic 
                mat_shader_name[0] = '\0';
            }

            char *Standard_Mattes = strstr(mat_name, ".Standard_Mattes");
            if (Standard_Mattes != NULL)
                //   Sources.Materials.myLibrary_ref_library.myMaterialName <<chop>> .Standard_Mattes
                Standard_Mattes[0] = '\0';

            const char * prefix = "Sources.Materials.";
            const char *mat_prefix_seperator = strstr(mat_name, prefix);
            if (mat_prefix_seperator != NULL) {
                //   Sources.Materials. <<SNIP>> myLibrary_ref_library.myMaterialName
                const char *mat_name_start = mat_prefix_seperator + strlen(prefix) ;
                strncpy(mat_name, mat_name_start, strlen(mat_name_start) + 1);
            }

            char *mdl_separator = strchr(mat_name, '.');
            if (strip_ns && mdl_separator != NULL) {
                //   myLibrary_ref_library. <<SNIP>> myMaterialName
                mdl_separator[0] = '\0';
                char *mat_name_start = mdl_separator + 1;
                strncpy(mat_name, mat_name_start, strlen(mat_name_start) + 1);
            } 
            // Leaving us with just the material name. 

        }

        // For maya, you get something simpler, like namespace:my_material_sg.

        char *ns_separator = strchr(mat_name, ':');
        if (strip_ns && ns_separator != NULL) {
            //    namespace: <<SNIP>> my_material_sg
            ns_separator[0] = '\0';
            char *mat_name_start = ns_separator + 1;
            strncpy(mat_name, mat_name_start, strlen(mat_name_start) + 1);
        } 

        strcpy(mat_name_out, mat_name);
    }
}

void get_clean_names(const char *mat_full_name, const char *obj_full_name, char *obj_name_out, char *mdl_name_out, char *mat_name_out, bool strip_mat_ns) {
    if (obj_full_name != NULL)
        get_clean_object_name  ( obj_full_name, obj_name_out, mdl_name_out);
    if (mat_full_name != NULL)
        get_clean_material_name(mat_full_name, mat_name_out, strip_mat_ns);
}

void override_name(char * name, const char * override_name) {
    if (override_name != NULL) {
        if (override_name[0] != NULL) {
            size_t override_length = strlen(override_name);
            strncpy(name, override_name, override_length );
            name[override_length] = NULL;
        }
    }
}

float hash_to_float(uint32_t hash) {
    uint32_t mantissa = hash & (( 1 << 23) - 1);
    uint32_t exponent = (hash >> 23) & ((1 << 8) - 1);
    exponent = std::max(exponent, (uint32_t) 1);
    exponent = std::min(exponent, (uint32_t) 254);
    exponent = exponent << 23;
    uint32_t sign = (hash >> 31);
    sign = sign << 31;
    uint32_t float_bits = sign | exponent | mantissa;
    float f;
    std::memcpy(&f, &float_bits, 4);
    return f;
}

void hash_name_rgb(char * name, AtColor* out_color) {
    // This puts the float ID into the red channel, and the human-readable
    // versions into the G and B channels. 
    uint32_t m3hash = 0;
    MurmurHash3_x86_32(name, (uint32_t) strlen(name), 0, &m3hash);
    out_color->r = hash_to_float(m3hash);
    out_color->g = ((float) ((m3hash << 8)) /  (float) UINT32_MAX);
    out_color->b = ((float) ((m3hash << 16)) / (float) UINT32_MAX);
}

const char* get_objects_shader_name(AtArray* shaders, AtUInt32 index = 0) {
    AtNode* shaderNode = static_cast<AtNode*>(AiArrayGetPtr( shaders, index));
    return AiNodeGetName(shaderNode);
}

void hash_object_rgb(AtShaderGlobals* sg, bool strip_mat_ns,
					 AtColor * mdl_hash_clr, AtColor * obj_hash_clr, AtColor * mat_hash_clr, 
					 const char * mdl_override, const char * obj_override,  const char * mat_override) {
    const char *mat_full_name = "";
    const char *obj_full_name = "";
    char obj_name[MAX_STRING_LENGTH];
    char mdl_name[MAX_STRING_LENGTH];
    char mat_name[MAX_STRING_LENGTH];

    AtNode* object = sg->Op;
    obj_full_name = AiNodeGetName(object);
    bool temp = false;

    AtArray * shaders = AiNodeGetArray(object, "shader");
    if (shaders->nelements == 1) {
        // Only one shader. No cluster materials or ambiguity. 
        mat_full_name = get_objects_shader_name(shaders, 0);
    } else {
        // There is a cluster material. We hash the name of the shader that is currently evaluating. 
        mat_full_name = AiNodeGetName(sg->shader);
    }

    get_clean_names(mat_full_name, obj_full_name, obj_name, mdl_name, mat_name, strip_mat_ns);

    override_name(obj_name, obj_override);
    override_name(mdl_name, mdl_override);
    override_name(mat_name, mat_override);

    hash_name_rgb( mdl_name, mdl_hash_clr );
    hash_name_rgb( obj_name, obj_hash_clr );
    hash_name_rgb( mat_name, mat_hash_clr );
}




///////////////////////////////////////////////
//
//      CRYPTOMATTE UTILITIES
//
///////////////////////////////////////////////



void write_array_of_AOVs(AtShaderGlobals * sg, AtArray * names, AtColor *color) {
    for (AtUInt32 i=0; i < names->nelements; i++) {
        const char * aovName = AiArrayGetStr( names, i);
        if (aovName == NULL) {
            return;
        }
        AtRGBA aovColor = AiRGBtoRGBA( *color );
        aovColor.b = AiColorToGrey(sg->out_opacity);

        if (strlen(aovName) > 0) {
            AiAOVSetRGBA(sg, aovName, aovColor);
        } else {
            return;
        }
    }
}



void override_settings_with_controller(int default_depth, int default_point_verb, int* depth_out, bool* strip_mat_ns_out) {
    AtNode * renderOptions = AiUniverseGetOptions();
    AtNode * automatte_controller = NULL;

    bool strip_mat_ns = true;
    int cryptomatte_depth = default_depth; 
    int point_verb = default_point_verb; 

    void * bg_node_pointer =  AiNodeGetPtr(renderOptions, "background");
    automatte_controller = static_cast<AtNode*>(bg_node_pointer);

    if (AiNodeIs(automatte_controller, "uBasic_controller")){
        *strip_mat_ns_out = AiNodeGetBool(automatte_controller, "strip_material_namespaces");
        bool override_cryptomatte = AiNodeGetBool(automatte_controller, "override_cryptomatte");
        if (override_cryptomatte) {
            cryptomatte_depth = AiNodeGetInt(automatte_controller, "cryptomatte_depth");
            cryptomatte_depth = std::min(cryptomatte_depth, MAX_CRYPTOMATTE_DEPTH);
            cryptomatte_depth = std::max(cryptomatte_depth, 1); // No fewer than 1, please.

            point_verb = AiNodeGetInt(automatte_controller, "pointcloud_instance_verbosity");
            point_verb = std::min(point_verb, 2);
            point_verb = std::max(point_verb, 0);
            g_pointcloud_instance_verbosity = point_verb;
        }
    }

    if (depth_out == NULL)
        return;
    
    // round up. 
    if ( cryptomatte_depth % 2 == 0 ) {
        *depth_out = cryptomatte_depth/2;
    } else {
        *depth_out = (cryptomatte_depth + 1)/2;
    }
}





///////////////////////////////////////////////
//
//      Metadata Writing
//
///////////////////////////////////////////////


typedef std::map<std::string,float>             md_map_type ;
typedef std::map<std::string,float>::iterator   md_map_iterator_type;


void write_metadata(AtNode * driver, std::string cryptomatte_name, md_map_type * map) {   
    AtArray * orig_md = AiNodeGetArray( driver, "custom_attributes");
    const AtUInt32 orig_num_entries = orig_md == NULL ? 0 : orig_md->nelements;

    const AtUInt32 new_num_entries = 4; // number of new entries
    std::string metadata_hash, metadata_conv, metadata_name, metadata_manf; // the new entries
    
    AtArray * combined_md = AiArrayAllocate(orig_num_entries + new_num_entries, 1, AI_TYPE_STRING); //Does not need destruction

    std::string prefix("STRING cryptomatte/");
    bool prefix_legal = false; // prefix is legal if it's not already used. For instance, "STRING cryptomatte/0/"
    int prefix_num = 0;
    char prefix_num_c[32];
    while (!prefix_legal) {
        sprintf(prefix_num_c, "%d", prefix_num);
        std::string prefix_candidate = prefix + std::string(prefix_num_c) + std::string("/");
        prefix_num++;
        prefix_legal = true; // assume innocence,
        for (AtUInt32 i=0; i<orig_num_entries; i++) {
            std::string entry(AiArrayGetStr(orig_md, i));
            prefix_legal = entry.compare(0, prefix_candidate.size(), prefix_candidate) != 0;
            if (!prefix_legal) // prove guilt. 
                break; 
        }
        if (prefix_legal)
            prefix = prefix_candidate;
    }

    metadata_hash = prefix + std::string("hash MurmurHash3_32");
    metadata_conv = prefix + std::string("conversion uint32_to_float32");
    metadata_name = prefix + std::string("name ") + cryptomatte_name;
    metadata_manf = prefix + std::string("manifest ");

    md_map_iterator_type map_it = map->begin();
    const size_t map_entries = map->size();
    const size_t max_entries = 100000;
    size_t metadata_entries = map_entries;
    if (map_entries > max_entries) {
        AiMsgWarning("Cryptomatte: %d entries in manifest, limiting to %d", map_entries, max_entries);
        metadata_entries = max_entries;
    }

    metadata_manf.append("{");
    for (AtUInt32 i=0; i<metadata_entries; i++) {
        const char * name = map_it->first.c_str();
        float hash_value = map_it->second;
        ++map_it;

        uint32_t float_bits;
        std::memcpy(&float_bits, &hash_value, 4);
        char hex_chars[9];
        sprintf(hex_chars, "%08x", float_bits);

        std::string pair;
        pair.append("\"");
        pair.append(name);
        pair.append("\":\"");
        pair.append(hex_chars);
        pair.append("\"");
        if (i < map_entries-1)
            pair.append(",");
        metadata_manf.append(pair);
    }
    metadata_manf.append("}");
    
    for (AtUInt32 i=0; i<orig_num_entries; i++) {
        AiArraySetStr(combined_md, i, AiArrayGetStr(orig_md, i));
    }
    AiArraySetStr(combined_md, orig_num_entries + 0, metadata_manf.c_str());
    AiArraySetStr(combined_md, orig_num_entries + 1, metadata_hash.c_str());
    AiArraySetStr(combined_md, orig_num_entries + 2, metadata_conv.c_str());
    AiArraySetStr(combined_md, orig_num_entries + 3, metadata_name.c_str());

    AiNodeSetArray( driver, "custom_attributes", combined_md);
}

#define METADATA_FLAG "already_has_crypto_metadata"

bool metadata_needed(AtNode* driver) {
    if (driver != NULL) {
        if (AiNodeLookUpUserParameter(driver, METADATA_FLAG) == NULL)
            return true;
    }
    return false;
}

void metadata_set_unneeded(AtNode* driver) {
    if (driver == NULL)
        return;
    if (AiNodeLookUpUserParameter(driver, METADATA_FLAG) == NULL)
        AiNodeDeclare(driver, METADATA_FLAG, "constant BOOL");
}


void add_hash_to_map(char * obj_name, md_map_type * md_map) {
    if (obj_name == NULL)
        return;
    if (obj_name[0] == NULL)
        return;
    AtColor hash;
    std::string name_string = obj_name;
    if (md_map->count(name_string) == 0) {
        hash_name_rgb(obj_name, &hash);
        (*md_map)[name_string] = hash.r;
    }
}


void add_hash_to_map_const(const char * name_in, md_map_type * md_map) {
    if (name_in == NULL)
        return;
    if (name_in[0] == NULL)
        return;
    char obj_name[MAX_STRING_LENGTH];               
    strcpy(obj_name, name_in);
    add_hash_to_map(obj_name, md_map);
}


void build_and_write_metadata(AtNode* driver_asset, AtNode* driver_object, AtNode* driver_material, 
                              std::string aov_asset, std::string aov_object, std::string aov_material, 
                              bool strip_mat_ns) {
    const clock_t metadata_start_time = clock();

    const bool do_md_asset = metadata_needed(driver_asset);
    const bool do_md_object = metadata_needed(driver_object);
    const bool do_md_material = metadata_needed(driver_material);

    metadata_set_unneeded(driver_asset);
    metadata_set_unneeded(driver_object);
    metadata_set_unneeded(driver_material);

    if (!do_md_asset && !do_md_object && !do_md_material)
        return;

    md_map_type map_md_asset;
    md_map_type map_md_object;
    md_map_type map_md_material;

    AtNodeIterator * shape_iterator = AiUniverseGetNodeIterator(AI_NODE_SHAPE);
    while (!AiNodeIteratorFinished(shape_iterator)) {
        AtNode *node = AiNodeIteratorGetNext(shape_iterator);
        const char *obj_full_name = AiNodeGetName(node);

        char obj_name[MAX_STRING_LENGTH];
        char mdl_name[MAX_STRING_LENGTH];
        char mat_name[MAX_STRING_LENGTH];
        get_clean_object_name(obj_full_name, obj_name, mdl_name);

        if (do_md_asset)
            add_hash_to_map(mdl_name, &map_md_asset);
        if (do_md_object)
            add_hash_to_map(obj_name, &map_md_object);

        if (do_md_material) {
            // Add all shaders from the object to the manifest. This includes cluster materials.
            AtArray * shaders = AiNodeGetArray(node, "shader");
            for (AtUInt32 i = 0; i < shaders->nelements; i++) {
                get_clean_material_name(get_objects_shader_name(shaders, i), mat_name, strip_mat_ns);
                add_hash_to_map(mat_name, &map_md_material); 
            }
        }
    }
    AiNodeIteratorDestroy(shape_iterator);

    AtNodeIterator * shader_iterator = AiUniverseGetNodeIterator(AI_NODE_SHADER);
    while (!AiNodeIteratorFinished(shader_iterator)) {
        AtNode *node = AiNodeIteratorGetNext(shader_iterator);
        if (AiNodeIs(node, "uBasic") || AiNodeIs(node, "automatte")) {
            const char *mat_full_name = AiNodeGetName(node);

            const char * override_asset = AiNodeGetStr( node, "override_asset");
            const char * override_object = AiNodeGetStr( node, "override_object");
            const char * override_material = AiNodeGetStr( node, "override_material");

            add_hash_to_map_const(override_asset, &map_md_asset);
            add_hash_to_map_const(override_object, &map_md_object);
            add_hash_to_map_const(override_material, &map_md_material);
        }
    }
    AiNodeIteratorDestroy(shader_iterator);

    write_metadata(driver_asset, aov_asset, &map_md_asset);
    write_metadata(driver_object, aov_object, &map_md_object);
    write_metadata(driver_material, aov_material, &map_md_material);

    AiMsgInfo("Cryptomatte manifest created - %f seconds", (float( clock () - metadata_start_time ) /  CLOCKS_PER_SEC));
}



///////////////////////////////////////////////
//
//      Building Cryptomatte Arnold Nodes
//
///////////////////////////////////////////////


void create_cryptomatte_aovs_filters_and_outputs(AutomatteData * data) {
    AtNode * renderOptions = AiUniverseGetOptions();

    AtArray * outputs = AiNodeGetArray( renderOptions, "outputs");

    int cryptomatte_aov_depth;
    bool strip_mat_ns = true;
    override_settings_with_controller(6, 0, &cryptomatte_aov_depth, &strip_mat_ns);

    AtArray * new_outputs = AiArrayAllocate(cryptomatte_aov_depth * 3, 1, AI_TYPE_STRING); // destroyed later
    int new_output_num = 0;

    AtNode * driver_cryptoAsset = NULL;
    AtNode * driver_cryptoObject = NULL;
    AtNode * driver_cryptoMaterial = NULL;

    for (AtUInt32 i=0; i < outputs->nelements; i++) {
        const char * output_string = AiArrayGetStr( outputs, i);

        size_t output_string_chars = strlen(output_string);

        char temp_string[MAX_STRING_LENGTH * 8]; 
        memset(temp_string, 0, sizeof(temp_string));
        strncpy(temp_string, output_string, output_string_chars);

        char * aov_name;
        char * aov_type_name;
        char * filter_name;
        char * driver_name;

        aov_name = strtok (temp_string," "); 
        aov_type_name = strtok (NULL," "); 
        filter_name = strtok (NULL," "); 
        driver_name = strtok (NULL," "); 

        size_t aov_len = strlen(aov_name);

        size_t aov_cryptoAsset_len = strlen(data->aov_cryptoasset.c_str());
        size_t aov_cryptoObject_len = strlen(data->aov_cryptoobject.c_str());
        size_t aov_cryptoMaterial_len = strlen(data->aov_cryptomaterial.c_str());

        AtNode * driver = NULL;
        AtArray * cryptoAOVs = NULL;
        if (aov_len == aov_cryptoAsset_len) {
            if (strncmp( aov_name, data->aov_cryptoasset.c_str(), aov_len) == 0) {
                cryptoAOVs = data->aovArray_cryptoasset;
                driver = AiNodeLookUpByName(driver_name);
                driver_cryptoAsset = driver;
            }
        }
        if (aov_len == aov_cryptoObject_len) {
            if (strncmp( aov_name, data->aov_cryptoobject.c_str(), aov_len) == 0) {
                cryptoAOVs = data->aovArray_cryptoobject;
                driver = AiNodeLookUpByName(driver_name);
                driver_cryptoObject = driver;

            }
        }
        if (aov_len == aov_cryptoMaterial_len) {
            if (strncmp( aov_name, data->aov_cryptomaterial.c_str(), aov_len) == 0) {
                cryptoAOVs = data->aovArray_cryptoMaterial;
                driver = AiNodeLookUpByName(driver_name);
                driver_cryptoMaterial = driver;

            }
        }


        if (cryptoAOVs != NULL) {

            ///////////////////////////////////////////////
            //
            //      Safety checks 
            //
            ///////////////////////////////////////////////

            if (!AiNodeIs(driver, "driver_exr")) {
                AiMsgError("Cryptomatte Error: Can only write Cryptomatte to EXR files.");
                return;
            }

            ///////////////////////////////////////////////
            //
            //      Compile info about original filter 
            //
            ///////////////////////////////////////////////

            float aFilter_width = 2.0;
            char aFilter_filter[128];
            AtNode * orig_filter = AiNodeLookUpByName(filter_name);
            const AtNodeEntry * orig_filter_nodeEntry = AiNodeGetNodeEntry(orig_filter);
            const char * orig_filter_type_name = AiNodeEntryGetName(orig_filter_nodeEntry);
            if (AiNodeEntryLookUpParameter(orig_filter_nodeEntry, "width") != NULL) {
                aFilter_width = AiNodeGetFlt(orig_filter, "width");             
            }

            memset(aFilter_filter, 0, sizeof(aFilter_filter));
            size_t filter_name_len = strlen(orig_filter_type_name);
            strncpy(aFilter_filter, orig_filter_type_name, filter_name_len);            
            char *filter_strip_point = strstr(aFilter_filter, "_filter");
            if (filter_strip_point != NULL) {
                filter_strip_point[0] = '\0';
            }

            ///////////////////////////////////////////////
            //
            //      Set CryptoAOV driver to full precision and outlaw RLE
            //
            ///////////////////////////////////////////////

            AiNodeSetBool(driver, "half_precision", false);

            const AtNodeEntry* driver_entry = AiNodeGetNodeEntry(driver);
            AtEnum compressions =  AiParamGetEnum(AiNodeEntryLookUpParameter(driver_entry, "compression"));         
            if (AiNodeGetInt(driver, "compression") == AiEnumGetValue(compressions, "rle")) {
                AiMsgWarning("Cryptomatte cannot be set to RLE compression- it does not work on full float. Switching to Zip.");
                AiNodeSetStr(driver, "compression", "zip");
            }
            
            ///////////////////////////////////////////////
            //
            //      Create filters and outputs as needed 
            //
            ///////////////////////////////////////////////

            for (int g=0; g<cryptomatte_aov_depth; g++) {
                char filter_rank_name[MAX_STRING_LENGTH];
                memset(filter_rank_name, 0, sizeof(filter_rank_name));

                char aov_rank_name[MAX_STRING_LENGTH];
                memset(aov_rank_name, 0, sizeof(aov_rank_name));

                size_t aov_name_chars = strlen(aov_name);
                strncpy(filter_rank_name, aov_name, aov_name_chars);
                strncpy(aov_rank_name, aov_name, aov_name_chars);

                char rank_number_string[MAX_STRING_LENGTH];
                memset(rank_number_string, 0, sizeof(rank_number_string));
                sprintf(rank_number_string, "%002d", g);

                strcat(filter_rank_name, "_filter" );
                strcat(filter_rank_name, rank_number_string );
                strcat(aov_rank_name, rank_number_string );
                
                if ( AiNodeLookUpByName( filter_rank_name ) == NULL) {
                    AtNode *filter = AiNode("cryptomatte_filter");
                    AiNodeSetStr(filter, "name", filter_rank_name);
                    AiNodeSetInt(filter, "rank", g*2);
                    AiNodeSetStr(filter, "filter", aFilter_filter);
                    AiNodeSetFlt(filter, "width", aFilter_width);
                    AiNodeSetStr(filter, "mode", "double_rgba");

                    AiAOVRegister(aov_rank_name, AI_TYPE_RGB, AI_AOV_BLEND_NONE);

                    // Add an output to the render globals, or make a list of outputs to add, and register an AOV
                    char new_output_string[MAX_STRING_LENGTH * 8];
                    memset(new_output_string, 0, sizeof(new_output_string));

                    strcat(new_output_string, aov_rank_name );
                    strcat(new_output_string, " " );
                    strcat(new_output_string, "RGBA" );
                    strcat(new_output_string, " " );
                    strcat(new_output_string, filter_rank_name );
                    strcat(new_output_string, " " );
                    strcat(new_output_string, driver_name );

                    AiArraySetStr(new_outputs, new_output_num, new_output_string);
                    new_output_num++;
                }

                AiArraySetStr(cryptoAOVs, g, aov_rank_name);
            }
        }
    }

    if (new_output_num > 0) {
        int total_outputs = outputs->nelements + new_output_num;
        AtArray * final_outputs = AiArrayAllocate(total_outputs, 1, AI_TYPE_STRING); // Does not need destruction
        for (AtUInt32 i=0; i < outputs->nelements; i++) {
            // Iterate through old outputs and add them
            AiArraySetStr(final_outputs, i, AiArrayGetStr( outputs, i));
        }
        for (int i=0; i < new_output_num; i++)  {
            // Iterate through new outputs and add them
            AiArraySetStr(final_outputs, i + outputs->nelements, AiArrayGetStr( new_outputs, i));
        }
        AiNodeSetArray(renderOptions, "outputs", final_outputs );
    }

    AiArrayDestroy(new_outputs);

    build_and_write_metadata(driver_cryptoAsset, driver_cryptoObject, driver_cryptoMaterial,
                             data->aov_cryptoasset, data->aov_cryptoobject, data->aov_cryptomaterial,
                             strip_mat_ns);
}


void create_cryptomatte_aovs_filters_and_outputs(const std::string& aov_cryptoasset, const std::string& aov_cryptoobject, const std::string& aov_cryptomaterial,
                                                   AtArray* aovArray_cryptoasset, AtArray* aovArray_cryptoobject, AtArray* aovArray_cryptoMaterial) {
    AtNode * renderOptions = AiUniverseGetOptions();

    AtArray * outputs = AiNodeGetArray( renderOptions, "outputs");

    int cryptomatte_aov_depth;
    bool strip_mat_ns = true;
    override_settings_with_controller(6, 0, &cryptomatte_aov_depth, &strip_mat_ns);

    AtArray * new_outputs = AiArrayAllocate(cryptomatte_aov_depth * 3, 1, AI_TYPE_STRING); // destroyed later
    int new_output_num = 0;

    AtNode * driver_cryptoAsset = NULL;
    AtNode * driver_cryptoObject = NULL;
    AtNode * driver_cryptoMaterial = NULL;

    for (AtUInt32 i=0; i < outputs->nelements; i++) {
        const char * output_string = AiArrayGetStr( outputs, i);

        size_t output_string_chars = strlen(output_string);

        char temp_string[MAX_STRING_LENGTH * 8]; 
        memset(temp_string, 0, sizeof(temp_string));
        strncpy(temp_string, output_string, output_string_chars);

        char * aov_name;
        char * aov_type_name;
        char * filter_name;
        char * driver_name;

        aov_name = strtok (temp_string," "); 
        aov_type_name = strtok (NULL," "); 
        filter_name = strtok (NULL," "); 
        driver_name = strtok (NULL," "); 

        size_t aov_len = strlen(aov_name);

        size_t aov_cryptoAsset_len = strlen(aov_cryptoasset.c_str());
        size_t aov_cryptoObject_len = strlen(aov_cryptoobject.c_str());
        size_t aov_cryptoMaterial_len = strlen(aov_cryptomaterial.c_str());

        AtNode * driver = NULL;
        AtArray * cryptoAOVs = NULL;
        if (aov_len == aov_cryptoAsset_len) {
            if (strncmp( aov_name, aov_cryptoasset.c_str(), aov_len) == 0) {
                cryptoAOVs = aovArray_cryptoasset;
                driver = AiNodeLookUpByName(driver_name);
                driver_cryptoAsset = driver;
            }
        }
        if (aov_len == aov_cryptoObject_len) {
            if (strncmp( aov_name, aov_cryptoobject.c_str(), aov_len) == 0) {
                cryptoAOVs = aovArray_cryptoobject;
                driver = AiNodeLookUpByName(driver_name);
                driver_cryptoObject = driver;

            }
        }
        if (aov_len == aov_cryptoMaterial_len) {
            if (strncmp( aov_name, aov_cryptomaterial.c_str(), aov_len) == 0) {
                cryptoAOVs = aovArray_cryptoMaterial;
                driver = AiNodeLookUpByName(driver_name);
                driver_cryptoMaterial = driver;

            }
        }


        if (cryptoAOVs != NULL) {

            ///////////////////////////////////////////////
            //
            //      Safety checks 
            //
            ///////////////////////////////////////////////

            if (!AiNodeIs(driver, "driver_exr")) {
                AiMsgError("Cryptomatte Error: Can only write Cryptomatte to EXR files.");
                return;
            }

            ///////////////////////////////////////////////
            //
            //      Compile info about original filter 
            //
            ///////////////////////////////////////////////

            float aFilter_width = 2.0;
            char aFilter_filter[128];
            AtNode * orig_filter = AiNodeLookUpByName(filter_name);
            const AtNodeEntry * orig_filter_nodeEntry = AiNodeGetNodeEntry(orig_filter);
            const char * orig_filter_type_name = AiNodeEntryGetName(orig_filter_nodeEntry);
            if (AiNodeEntryLookUpParameter(orig_filter_nodeEntry, "width") != NULL) {
                aFilter_width = AiNodeGetFlt(orig_filter, "width");             
            }

            memset(aFilter_filter, 0, sizeof(aFilter_filter));
            size_t filter_name_len = strlen(orig_filter_type_name);
            strncpy(aFilter_filter, orig_filter_type_name, filter_name_len);            
            char *filter_strip_point = strstr(aFilter_filter, "_filter");
            if (filter_strip_point != NULL) {
                filter_strip_point[0] = '\0';
            }

            ///////////////////////////////////////////////
            //
            //      Set CryptoAOV driver to full precision and outlaw RLE
            //
            ///////////////////////////////////////////////

            AiNodeSetBool(driver, "half_precision", false);

            const AtNodeEntry* driver_entry = AiNodeGetNodeEntry(driver);
            AtEnum compressions =  AiParamGetEnum(AiNodeEntryLookUpParameter(driver_entry, "compression"));         
            if (AiNodeGetInt(driver, "compression") == AiEnumGetValue(compressions, "rle")) {
                AiMsgWarning("Cryptomatte cannot be set to RLE compression- it does not work on full float. Switching to Zip.");
                AiNodeSetStr(driver, "compression", "zip");
            }
            
            ///////////////////////////////////////////////
            //
            //      Create filters and outputs as needed 
            //
            ///////////////////////////////////////////////

            for (int g=0; g<cryptomatte_aov_depth; g++) {
                char filter_rank_name[MAX_STRING_LENGTH];
                memset(filter_rank_name, 0, sizeof(filter_rank_name));

                char aov_rank_name[MAX_STRING_LENGTH];
                memset(aov_rank_name, 0, sizeof(aov_rank_name));

                size_t aov_name_chars = strlen(aov_name);
                strncpy(filter_rank_name, aov_name, aov_name_chars);
                strncpy(aov_rank_name, aov_name, aov_name_chars);

                char rank_number_string[MAX_STRING_LENGTH];
                memset(rank_number_string, 0, sizeof(rank_number_string));
                sprintf(rank_number_string, "%002d", g);

                strcat(filter_rank_name, "_filter" );
                strcat(filter_rank_name, rank_number_string );
                strcat(aov_rank_name, rank_number_string );
                
                if ( AiNodeLookUpByName( filter_rank_name ) == NULL) {
                    AtNode *filter = AiNode("cryptomatte_filter");
                    AiNodeSetStr(filter, "name", filter_rank_name);
                    AiNodeSetInt(filter, "rank", g*2);
                    AiNodeSetStr(filter, "filter", aFilter_filter);
                    AiNodeSetFlt(filter, "width", aFilter_width);
                    AiNodeSetStr(filter, "mode", "double_rgba");

                    AiAOVRegister(aov_rank_name, AI_TYPE_RGB, AI_AOV_BLEND_NONE);

                    // Add an output to the render globals, or make a list of outputs to add, and register an AOV
                    char new_output_string[MAX_STRING_LENGTH * 8];
                    memset(new_output_string, 0, sizeof(new_output_string));

                    strcat(new_output_string, aov_rank_name );
                    strcat(new_output_string, " " );
                    strcat(new_output_string, "RGBA" );
                    strcat(new_output_string, " " );
                    strcat(new_output_string, filter_rank_name );
                    strcat(new_output_string, " " );
                    strcat(new_output_string, driver_name );

                    AiArraySetStr(new_outputs, new_output_num, new_output_string);
                    new_output_num++;
                }

                AiArraySetStr(cryptoAOVs, g, aov_rank_name);
            }
        }
    }

    if (new_output_num > 0) {
        int total_outputs = outputs->nelements + new_output_num;
        AtArray * final_outputs = AiArrayAllocate(total_outputs, 1, AI_TYPE_STRING); // Does not need destruction
        for (AtUInt32 i=0; i < outputs->nelements; i++) {
            // Iterate through old outputs and add them
            AiArraySetStr(final_outputs, i, AiArrayGetStr( outputs, i));
        }
        for (int i=0; i < new_output_num; i++)  {
            // Iterate through new outputs and add them
            AiArraySetStr(final_outputs, i + outputs->nelements, AiArrayGetStr( new_outputs, i));
        }
        AiNodeSetArray(renderOptions, "outputs", final_outputs );
    }

    AiArrayDestroy(new_outputs);

    build_and_write_metadata(driver_cryptoAsset, driver_cryptoObject, driver_cryptoMaterial,
                             aov_cryptoasset, aov_cryptoobject, aov_cryptomaterial,
                             strip_mat_ns);
}


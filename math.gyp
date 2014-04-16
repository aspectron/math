{
    'variables': {
        # should point to jsx folder. Also used in the included gypi files
        'JSX' : '../jsx-priv/',
    }, # variables
    'includes': [
        'build/common.gypi',
        'build/debug_x86.gypi',
        'build/debug_x64.gypi',
        'build/release_x86.gypi',
        'build/release_x64.gypi',
    ], # includes

    'target_defaults': {
        'msvs_settings': {
            'VCCLCompilerTool': {
                'PreprocessorDefinitions': [
                    'MATH_EXPORTS',
                ], # PreprocessorDefinitions
            }, # VCCLCompilerTool
        }, # msvs_settings

        'configurations': {
            'Debug_Win32': {
                'inherit_from' : ['Debug_x86'],
            }, # Debug_Win32

            'Debug_x64': {
                'inherit_from' : ['Debug64'],
            }, # Debug_x64

            'Release_Win32': {
                'inherit_from' : ['Release_x86'],
            }, # Release_Win32

            'Release_x64': {
                'inherit_from' : ['Release64'],
            }, # Release_x64
        }, # configurations

        'sources': [
            'src/math.euler_angles.cpp',
            'src/math.matrix.cpp',
            'src/math.noise.cpp',
            'src/math.orientation_vectors.cpp',
            'src/math.quaternion.cpp',
            'src/math.vector.cpp',
            'src/math.euler_angles.hpp',
            'src/math.hpp',
            'src/math.matrix.hpp',
            'src/math.noise.hpp',
            'src/math.orientation_vectors.hpp',
            'src/math.quaternion.hpp',
            'src/math.vector.hpp',
            'src/math.vector2.hpp',
            'src/math.vector3.hpp',
            'src/math.vector4.hpp',
        ], # sources
    }, # target_defaults

    'targets': [
        {
            'target_name': 'math',
            'type': 'shared_library',
        }, # target - math
    ], # targets
}
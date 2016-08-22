{
    'variables': {
        'include_files': [
            'include/math/math.hpp',
            'include/math/euler_angles.hpp',
            'include/math/matrix.hpp',
            'include/math/noise.hpp',
            'include/math/orientation_vectors.hpp',
            'include/math/quaternion.hpp',
            'include/math/vector2.hpp',
            'include/math/vector3.hpp',
            'include/math/vector4.hpp',
        ],
        'source_files': [
            'src/math.cpp',
            'src/euler_angles.cpp',
            'src/matrix.cpp',
            'src/noise.cpp',
            'src/orientation_vectors.cpp',
            'src/quaternion.cpp',
        ],
        'include_dirs': ['include', '<!(node -e require(\'v8pp\'))'],
    },
    'targets': [
        {
            'target_name': 'math',
            'cflags_cc+': ['-std=c++11', '-fexceptions'],
            'msvs_settings': { 'VCCLCompilerTool': { 'ExceptionHandling': 1 } },
            'xcode_settings': { 'GCC_ENABLE_CPP_EXCEPTIONS': 'YES' },
            'include_dirs': ['<@(include_dirs)'],
            'direct_dependent_settings': {
                'include_dirs': ['<@(include_dirs)'],
            },
            'defines': ['MATH_EXPORTS'],
            'defines!': ['V8_DEPRECATION_WARNINGS=1'],
            'sources': ['<@(include_files)', '<@(source_files)'],
        },
    ],
}
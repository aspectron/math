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
            'src/euler_angles.cpp',
            'src/matrix.cpp',
            'src/noise.cpp',
            'src/orientation_vectors.cpp',
            'src/quaternion.cpp',
        ],
    },
    'targets': [
        {
            'target_name': 'math',
            'type': 'shared_library',
            'msvs_guid': 'AF127E22-2727-482C-9C41-4F7575F984DC',
            'dependencies': [
                '<(jsx)/jsx-lib.gyp:jsx-lib',
                '<(jsx)/extern/extern.gyp:*',
            ],
            'include_dirs': ['include'],
            'direct_dependent_settings': {
                'include_dirs': ['include'],
            },
            'defines': ['MATH_EXPORTS'],
            'sources': ['<@(include_files)', '<@(source_files)'],
        },
    ],
}
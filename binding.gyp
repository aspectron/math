{
    'includes': ['common.gypi'],
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
            'include_dirs': ['<@(include_dirs)'],
            'direct_dependent_settings': {
                'include_dirs': ['<@(include_dirs)'],
            },
            'defines': ['MATH_EXPORTS'],
            'sources': ['<@(include_files)', '<@(source_files)'],
        },
    ],
}
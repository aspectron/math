{
    'targets': [
        {
            'target_name': 'math',
            'type': 'shared_library',
            'msvs_guid': 'AF127E22-2727-482C-9C41-4F7575F984DC',
            'dependencies': [
                '<(jsx)/sdk/core/core.gyp:core',
                '<(jsx)/extern/extern.gyp:*',
            ],
            'direct_dependent_settings': {
                'include_dirs': ['src'],
            },
            'defines': ['MATH_EXPORTS'],
            'sources': [
                'src/math.hpp',
                'src/math.euler_angles.cpp',
                'src/math.euler_angles.hpp',
                'src/math.matrix.cpp',
                'src/math.matrix.hpp',
                'src/math.noise.cpp',
                'src/math.noise.hpp',
                'src/math.orientation_vectors.cpp',
                'src/math.orientation_vectors.hpp',
                'src/math.quaternion.cpp',
                'src/math.quaternion.hpp',
                'src/math.vector.cpp',
                'src/math.vector.hpp',
                'src/math.vector2.hpp',
                'src/math.vector3.hpp',
                'src/math.vector4.hpp',
            ],
        },
    ],
}
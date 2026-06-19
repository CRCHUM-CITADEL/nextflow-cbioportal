process BUILD_ANON_LINKING {
    label 'python'

    input:
    path(linking_file)

    output:
    path("anon_linking.txt")

    script:
    """
    build_anon_linking.py ${linking_file}
    """
}

@PKGNAME@ <- local(get(load(system.file(file.path("extdata", "@VCFDBDATA@"), package="@PKGNAME@"))))
@PKGNAME@@db.file <- system.file(file.path("extdata", "@DB_NAME@"), package="@PKGNAME@")

#@PKGNAME@ <- new("VcfDB", db.path=system.file(file.path("extdata", "@DB_NAME@"), package="@PKGNAME@"), tbsl=tbsl)

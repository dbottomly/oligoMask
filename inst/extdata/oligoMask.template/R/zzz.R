tbsl <- local(get(load(system.file(file.path("extdata", "@TBSLDATA@"), package="@PKGNAME@"))))
@PKGNAME@ <- new("VcfDB", db.path=system.file(file.path("extdata", "@DB_NAME@"), package="@PKGNAME@"), tbsl=tbsl)

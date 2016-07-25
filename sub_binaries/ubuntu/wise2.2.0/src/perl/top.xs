
#ifdef __cplusplus
extern "C" {
#endif
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#ifdef __cplusplus
}
#endif

#include "wise2.h"


static int
not_here(s)
char *s;
{
    croak("%s not implemented on this architecture", s);
    return -1;
}

static double
constant(name, arg)
char *name;
int arg;
{
    errno = 0;
    switch (*name) {
    }
    errno = EINVAL;
    return 0;

not_there:
    errno = ENOENT;
    return 0;
}


MODULE = Wise2 PACKAGE = Wise2


double
constant(name,arg)
	char *		name
	int		arg


PROTOTYPES: ENABLE

void
error_on(type)
	int type;
	CODE:
	Wise2_error_on(type);


void
error_off(type)
	int type;
	CODE:
	Wise2_error_off(type);

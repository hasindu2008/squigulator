/* @file  misc.c
**
** @@
******************************************************************************/

#include "error.h"

enum sq_log_level_opt _log_level = LOG_VERB;

enum sq_log_level_opt get_log_level(){
    return _log_level;
}

void set_log_level(enum sq_log_level_opt level){
    _log_level = level;
}

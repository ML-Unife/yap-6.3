/*************************************************************************
*									 *
*	 YAP Prolog 							 *
*									 *
*	Yap Prolog was developed at NCCUP - Universidade do Porto	 *
*									 *
* Copyright L.Damas, V.S.Costa and Universidade do Porto 1985-1997	 *
*									 *
**************************************************************************
*									 *
* File:		load_dl.c						 *
* comments:	dl based dynamic loaderr of external routines		 *
*               tested on i486-linuxelf					 *
*************************************************************************/

#include "Yap.h"
#include "Yatom.h"
#include "Heap.h"
#include "yapio.h"
#include "Foreign.h"

#if LOAD_DL

#include <dlfcn.h>
#include <string.h>
#include <stdio.h>

/*
 *   YAP_FindExecutable(argv[0]) should be called on yap initialization to
 *   locate the executable of Yap
*/
void
Yap_FindExecutable(char *name)
{
}


/*
 * LoadForeign(ofiles,libs,proc_name,init_proc) dynamically loads foreign
 * code files and libraries and locates an initialization routine
*/
static Int
LoadForeign(StringList ofiles, StringList libs,
	       char *proc_name,	YapInitProc *init_proc)
{

  while (ofiles) {
    void *handle;

    /* dlopen wants to follow the LD_CONFIG_PATH */
    if (!Yap_TrueFileName(ofiles->s, Yap_FileNameBuf, TRUE)) {
      strcpy(Yap_ErrorSay, "[ Trying to open unexisting file in LoadForeign ]");
      return LOAD_FAILLED;
    }
#ifdef __osf__
    if((handle=dlopen(Yap_FileNameBuf,RTLD_LAZY)) == 0)
#else
    if((handle=dlopen(Yap_FileNameBuf,RTLD_LAZY|RTLD_GLOBAL)) == 0)
#endif
    {
      fprintf(stderr,"dlopen of %s failed with error %s\n", Yap_FileNameBuf, dlerror());
/*      strcpy(Yap_ErrorSay,dlerror());*/
      return LOAD_FAILLED;
    }

    ofiles->handle = handle;

    if (!*init_proc)
      *init_proc = (YapInitProc) dlsym(handle,proc_name);

    ofiles = ofiles->next;
  }

  if(! *init_proc) {
    strcpy(Yap_ErrorSay,"Could not locate initialization routine");
    return LOAD_FAILLED;
  }

  /* load libraries first so that their symbols are available to
     other routines */
  while (libs) {

    if (libs->s[0] == '-') {
      strcpy(Yap_FileNameBuf,"lib");
      strcat(Yap_FileNameBuf,libs->s+2);
      strcat(Yap_FileNameBuf,".so");
    } else {
      strcpy(Yap_FileNameBuf,libs->s);
    }

#ifdef __osf__
    if((libs->handle=dlopen(Yap_FileNameBuf,RTLD_LAZY)) == NULL)
#else
    if((libs->handle=dlopen(Yap_FileNameBuf,RTLD_LAZY|RTLD_GLOBAL)) == NULL)
#endif
    {
      strcpy(Yap_ErrorSay,dlerror());
      return LOAD_FAILLED;
    }
    libs = libs->next;
  }
  return LOAD_SUCCEEDED;
}

Int
Yap_LoadForeign(StringList ofiles, StringList libs,
	       char *proc_name,	YapInitProc *init_proc)
{
  return LoadForeign(ofiles, libs, proc_name, init_proc);
}

void 
Yap_ShutdownLoadForeign(void)
{
  ForeignObj *f_code;

  f_code = ForeignCodeLoaded;
  while (f_code != NULL) {
    StringList objs, libs;

    objs = f_code->objs;
    while (objs != NULL) {
      if (dlclose(objs->handle) != 0)
	return; /* ERROR */
      objs = objs->next;
    }
    libs = f_code->libs;
    while (libs != NULL) {
      if (dlclose(libs->handle) != 0)
	return; /* ERROR */
      objs = libs->next;
    }
    f_code = f_code->next;
  }
}

Int
Yap_ReLoadForeign(StringList ofiles, StringList libs,
	       char *proc_name,	YapInitProc *init_proc)
{
  return(LoadForeign(ofiles,libs, proc_name, init_proc));
}

#endif

#if SIMICS

void dlopen(void)
{
}

void dlclose(void)
{
}

void dlsym(void)
{
}

#endif




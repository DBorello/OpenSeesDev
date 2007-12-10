/* 
 * tclMain.c --
 *
 *	Main program for Tcl shells and other Tcl-based applications.
 *
 * Copyright (c) 1988-1994 The Regents of the University of California.
 * Copyright (c) 1994-1997 Sun Microsystems, Inc.
 *
 * See the file "license.terms" for information on usage and redistribution
 * of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 *
 * RCS: @(#) $Id: tclMain.cpp,v 1.43 2007-12-10 23:31:06 fmk Exp $
 */

/*                       MODIFIED   FOR                              */

/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

#include <string.h>

#ifndef _WIN32
#include <unistd.h>
#endif

extern "C" {
#include <tcl.h>
#include <tclDecls.h>
EXTERN int		TclFormatInt _ANSI_ARGS_((char *buffer, long n));
EXTERN int		TclObjCommandComplete _ANSI_ARGS_((Tcl_Obj *cmdPtr));
}


#ifdef _PARALLEL_INTERPRETERS
#include <MachineBroker.h>
extern MachineBroker *theMachineBroker;
#endif

const char * getInterpPWD(Tcl_Interp *interp);

#include <OPS_Globals.h>

int		Tcl_AppInit _ANSI_ARGS_((Tcl_Interp *interp));


# undef TCL_STORAGE_CLASS
# define TCL_STORAGE_CLASS DLLEXPORT

/*
 * The following code ensures that tclLink.c is linked whenever
 * Tcl is linked.  Without this code there's no reference to the
 * code in that file from anywhere in Tcl, so it may not be
 * linked into the application.
 */

#ifdef _TCL84
int (*tclDummyLinkVarPtr)(Tcl_Interp *interp, const char *a,
			  char *b, int c) = Tcl_LinkVar;
#else
int (*tclDummyLinkVarPtr)(Tcl_Interp *interp, char *a,
			  char *b, int c) = Tcl_LinkVar;
#endif

/*
 * Declarations for various library procedures and variables (don't want
 * to include tclPort.h here, because people might copy this file out of
 * the Tcl source directory to make their own modified versions).
 * Note:  "exit" should really be declared here, but there's no way to
 * declare it without causing conflicts with other definitions elsewher
 * on some systems, so it's better just to leave it out.
 */


typedef struct parameterValues {
  char *value;
  struct parameterValues *next;
} OpenSeesTcl_ParameterValues;

typedef struct parameter {
  char *name;
  OpenSeesTcl_ParameterValues *values;
  struct parameter *next;
} OpenSeesTcl_Parameter;


#ifdef _WIN32
extern "C" int	isatty _ANSI_ARGS_((int fd));
extern "C" char * strcpy _ANSI_ARGS_((char *dst, CONST char *src)) throw();
#endif
static char *tclStartupScriptFileName = NULL;


#include <FileStream.h>
#include <SimulationInformation.h>
SimulationInformation simulationInfo;

char *simulationInfoOutputFilename = 0;
char *neesCentralProjID =0;
char * neesCentralExpID =0;
char *neesCentralUser =0;
char *neesCentralPasswd =0;


/*
 *----------------------------------------------------------------------
 *
 * TclSetStartupScriptFileName --
 *
 *	Primes the startup script file name, used to override the
 *      command line processing.
 *
 * Results:
 *	None. 
 *
 * Side effects:
 *	This procedure initializes the file name of the Tcl script to
 *      run at startup.
 *
 *----------------------------------------------------------------------
 */
void TclSetStartupScriptFileName(char *fileName)
{
    tclStartupScriptFileName = fileName;
}


/*
 *----------------------------------------------------------------------
 *
 * TclGetStartupScriptFileName --
 *
 *	Gets the startup script file name, used to override the
 *      command line processing.
 *
 * Results:
 *	The startup script file name, NULL if none has been set.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */
char *TclGetStartupScriptFileName()
{
    return tclStartupScriptFileName;
}

int
EvalFileWithParameters(Tcl_Interp *interp, 
		       char *tclStartupFileScript, 
		       OpenSeesTcl_Parameter *theParameters, 
		       char **paramNames, 
		       char **paramValues, 
		       int numParam, 
		       int currentParam, 
		       int rank, 
		       int np)
{
  if (currentParam < numParam) {
    OpenSeesTcl_Parameter *theCurrentParam = theParameters;
    OpenSeesTcl_Parameter *theNextParam = theParameters->next;
    char *paramName = theCurrentParam->name;
    paramNames[currentParam] = paramName;

    OpenSeesTcl_ParameterValues *theValue = theCurrentParam->values;
    int nextParam = currentParam+1;
    while (theValue != 0) {
      char *paramValue = theValue->value;
      paramValues[currentParam] = paramValue;
      EvalFileWithParameters(interp, 
			     tclStartupFileScript, 
			     theNextParam, 
			     paramNames, 
			     paramValues, 
			     numParam, 
			     nextParam, 
			     rank, 
			     np);

      theValue=theValue->next;
    } 
  } else {
    
    simulationInfo.start();
    static int count = 0;
    
    if ((count % np) == rank) {
      Tcl_Eval(interp, "wipe");
     
	  
      for (int i=0; i<numParam; i++) {
		  
	Tcl_SetVar(interp, paramNames[i], paramValues[i], TCL_GLOBAL_ONLY);	    
	simulationInfo.addParameter(paramNames[i], paramValues[i]); 
     }

      count++;

      const char *pwd = getInterpPWD(interp);
      simulationInfo.addInputFile(tclStartupScriptFileName, pwd);

      int ok = Tcl_EvalFile(interp, tclStartupScriptFileName);

      simulationInfo.end();
      
      return ok;
    }
    else
      count++;
  }

  return 0;
}



/*
 *----------------------------------------------------------------------
 *
 * Tcl_Main --
 *
 *	Main program for tclsh and most other Tcl-based applications.
 *
 * Results:
 *	None. This procedure never returns (it exits the process when
 *	it's done.
 *
 * Side effects:
 *	This procedure initializes the Tcl world and then starts
 *	interpreting commands;  almost anything could happen, depending
 *	on the script being interpreted.
 *
 *----------------------------------------------------------------------
 */

//// Boris Jeremic {
//#include <fstream.h>
//// } Boris Jeremic


void
g3TclMain(int argc, char **argv, Tcl_AppInitProc * appInitProc, int rank, int np)
{ 
    Tcl_Obj *resultPtr;
    Tcl_Obj *commandPtr = NULL;
    char buffer[1000], *args;
    int code, gotPartial, tty, length;
    int exitCode = 0;
    Tcl_Channel inChannel, outChannel, errChannel;
    Tcl_Interp *interp;
    Tcl_DString argString;

#ifdef _PARALLEL_INTERPRETERS
    if (theMachineBroker->getPID() == 0) {
#endif

	/* fmk - beginning of modifications for OpenSees */
	fprintf(stderr,"\n\n\t OpenSees -- Open System For Earthquake Engineering Simulation");
	fprintf(stderr,"\n\tPacific Earthquake Engineering Research Center -- %s\n\n", OPS_VERSION);
	
	fprintf(stderr,"\t    (c) Copyright 1999,2000 The Regents of the University of California");
	fprintf(stderr,"\n\t\t\t\t All Rights Reserved\n");    
	fprintf(stderr,"    (Copyright and Disclaimer @ http://www.berkeley.edu/OpenSees/copyright.html)\n\n\n");
	

#ifdef _PARALLEL_INTERPRETERS
    }
#endif

    Tcl_FindExecutable(argv[0]);
    interp = Tcl_CreateInterp();


#ifdef TCL_MEM_DEBUG
    Tcl_InitMemory(interp);
#endif

    /*
     * Make command-line arguments available in the Tcl variables "argc"
     * and "argv".  If the first argument doesn't start with a "-" then
     * strip it off and use it as the name of a script file to process.
     */

    if (tclStartupScriptFileName == NULL) {
	if ((argc > 1) && (argv[1][0] != '-')) {
	    tclStartupScriptFileName = argv[1];
	    argc--;
	    argv++;
	}
    }

    args = Tcl_Merge(argc-1, argv+1);
    Tcl_ExternalToUtfDString(NULL, args, -1, &argString);
    Tcl_SetVar(interp, "argv", Tcl_DStringValue(&argString), TCL_GLOBAL_ONLY);
    Tcl_DStringFree(&argString);
    ckfree(args);


    if (tclStartupScriptFileName == NULL) {
	Tcl_ExternalToUtfDString(NULL, argv[0], -1, &argString);
    } else {
	tclStartupScriptFileName = Tcl_ExternalToUtfDString(NULL,
		tclStartupScriptFileName, -1, &argString);
    }

    TclFormatInt(buffer, argc-1);
    Tcl_SetVar(interp, "argc", buffer, TCL_GLOBAL_ONLY);
    Tcl_SetVar(interp, "argv0", Tcl_DStringValue(&argString), TCL_GLOBAL_ONLY);

    /*
     * Set the "tcl_interactive" variable.
     */

    tty = isatty(0);
    char one[2] = "1";
    char zero[2] = "0";

    Tcl_SetVar(interp, "tcl_interactive",
	    ((tclStartupScriptFileName == NULL) && tty) ? one : zero,
	    TCL_GLOBAL_ONLY);
    
    /*
     * Invoke application-specific initialization.
     */

    if ((*appInitProc)(interp) != TCL_OK) {
	errChannel = Tcl_GetStdChannel(TCL_STDERR);
	if (errChannel) {
	    Tcl_WriteChars(errChannel,
		    "application-specific initialization failed: ", -1);
	    Tcl_WriteObj(errChannel, Tcl_GetObjResult(interp));
	    Tcl_WriteChars(errChannel, "\n", 1);
	}
    }

    /*
     * If a script file was specified then just source that file
     * and quit.
     */


    if (tclStartupScriptFileName != NULL) {
      OpenSeesTcl_Parameter *theParameters = 0;
      OpenSeesTcl_Parameter *endParameters = 0;
      int numParam = 0;

      if (argc > 1) {
	int currentArg = 1;
	while (currentArg < argc && argv[currentArg] != NULL) {

	  if ((strcmp(argv[currentArg], "-par") == 0) || (strcmp(argv[currentArg], "-Par") == 0)) {
	    
	    if (argc > (currentArg+2)) {
	      
	      char *parName = argv[currentArg+1];
	      char *parValue = argv[currentArg+2];
	      
	      // add a OpenSeesTcl_Parameter to end of list of parameters
	      OpenSeesTcl_Parameter *nextParam = new OpenSeesTcl_Parameter;
	      nextParam->name = new char [strlen(parName)+1];
	      strcpy(nextParam->name, parName);
	      nextParam->values = 0;
	      
	      if (theParameters == 0)
		theParameters = nextParam;
	      if (endParameters != 0)
		endParameters->next = nextParam;
	      nextParam->next = 0;
	      endParameters = nextParam;
	      
	      // now open par values files to create the values
	      char nextLine[1000];
	      FILE *valueFP = fopen(parValue,"r");
	      if (valueFP != 0) {
		OpenSeesTcl_ParameterValues *endValues = 0;
		
		while (fscanf(valueFP, "%s", nextLine) != EOF) {
		  
		  OpenSeesTcl_ParameterValues *nextValue = new OpenSeesTcl_ParameterValues;
		  nextValue->value = new char [strlen(nextLine)+1];
		  strcpy(nextValue->value, nextLine);
		  
		  if (nextParam->values == 0) {
		    nextParam->values = nextValue;
		  }
		if (endValues != 0)
		  endValues->next = nextValue;
		endValues = nextValue;
		nextValue->next = 0;	      
		}
		fclose(valueFP);
	      } else {
		
		OpenSeesTcl_ParameterValues *nextValue = new OpenSeesTcl_ParameterValues;		
		nextValue->value = new char [strlen(parValue)+1];
		
		strcpy(nextValue->value, parValue);
		
		nextParam->values = nextValue;
		nextValue->next = 0;
		
	      }
	      numParam++;
	    }
	    currentArg += 3;
	  } else if ((strcmp(argv[currentArg], "-info") == 0) || (strcmp(argv[currentArg], "-INFO") == 0)) {
	    if (argc > (currentArg+1)) {
	      simulationInfoOutputFilename = argv[currentArg+1];	    
	    }			   
	    currentArg+=2;
	  } else if ((strcmp(argv[currentArg], "-upload") == 0) || (strcmp(argv[currentArg], "-UPLOAD") == 0)) {
	    bool more = true;
	    currentArg++;
	    while (more == true && currentArg < argc) {
	      
	      if (strcmp(argv[currentArg],"-user") == 0) {
		neesCentralUser = argv[currentArg+1];
		currentArg += 2;

	      } else if (strcmp(argv[currentArg],"-pass") == 0) {
		neesCentralPasswd = argv[currentArg+1];
		currentArg += 2;
	      } else if (strcmp(argv[currentArg],"-projID") == 0) {
		neesCentralProjID = argv[currentArg+1];
		currentArg += 2;
		
	      } else if (strcmp(argv[currentArg],"-expID") == 0) {
		neesCentralExpID = argv[currentArg+1];
		currentArg += 2;
    
	      } else if (strcmp(argv[currentArg],"-title") == 0) {
		simulationInfo.setTitle(argv[currentArg+1]);	
		currentArg += 2;
      
	      } else if (strcmp(argv[currentArg],"-description") == 0) {
		simulationInfo.setDescription(argv[currentArg+1]);	
		currentArg += 2;      
	      } else
		more = false;
	    }
	  } else 
	    currentArg++;
	}

	if (numParam != 0) {
	  char **paramNames = new char *[numParam];
	  char **paramValues = new char *[numParam];
	  code = EvalFileWithParameters(interp, tclStartupScriptFileName, theParameters, paramNames, paramValues, numParam, 0, rank, np);

	  if (code != TCL_OK) {
	    errChannel = Tcl_GetStdChannel(TCL_STDERR);
	    if (errChannel) {
	      /*
	       * The following statement guarantees that the errorInfo
	       * variable is set properly.
	       */
	      
	      Tcl_AddErrorInfo(interp, "");
	      Tcl_WriteObj(errChannel, Tcl_GetVar2Ex(interp, "errorInfo",
						     NULL, TCL_GLOBAL_ONLY));
	      Tcl_WriteChars(errChannel, "\n", 1);
	    }
	    exitCode = 1;
	  }
	  goto done;
	}
      }

      //      if (simulationInfoOutputFilename != 0) {
      const char *pwd = getInterpPWD(interp);
      simulationInfo.start();
      simulationInfo.addInputFile(tclStartupScriptFileName, pwd);
      
      code = Tcl_EvalFile(interp, tclStartupScriptFileName);
      
      if (code != TCL_OK) {
	errChannel = Tcl_GetStdChannel(TCL_STDERR);
	if (errChannel) {
	  /*
	   * The following statement guarantees that the errorInfo
	   * variable is set properly.
	   */
	  
	  Tcl_AddErrorInfo(interp, "");
	  Tcl_WriteObj(errChannel, Tcl_GetVar2Ex(interp, "errorInfo",
						 NULL, TCL_GLOBAL_ONLY));
	  Tcl_WriteChars(errChannel, "\n", 1);
	}
	exitCode = 1;
      }
      goto done;
    } else {

    /*
     * We're running interactively.  Source a user-specific startup
     * file if the application specified one and if the file exists.
     */

      Tcl_DStringFree(&argString);
      
      int currentArg = 1;
      while (currentArg < argc && argv[currentArg] != NULL) {
	if ((strcmp(argv[currentArg], "-info") == 0) || (strcmp(argv[currentArg], "-INFO") == 0)) {
	  if (argc > (currentArg+1)) {
	    simulationInfoOutputFilename = argv[currentArg+1];	    
	  }			   
	  currentArg+=2;
	} else if ((strcmp(argv[currentArg], "-upload") == 0) || (strcmp(argv[currentArg], "-UPLOAD") == 0)) {
	  bool more = true;
	  currentArg++;
	  while (more == true && currentArg < argc) {
	    
	    if (strcmp(argv[currentArg],"-user") == 0) {
	      neesCentralUser = argv[currentArg+1];
	      currentArg += 2;
	      
	    } else if (strcmp(argv[currentArg],"-pass") == 0) {
	      neesCentralPasswd = argv[currentArg+1];
	      currentArg += 2;
	    } else if (strcmp(argv[currentArg],"-projID") == 0) {
	      neesCentralProjID = argv[currentArg+1];
	      currentArg += 2;
	      
	    } else if (strcmp(argv[currentArg],"-expID") == 0) {
	      neesCentralExpID = argv[currentArg+1];
	      currentArg += 2;
	      
	    } else if (strcmp(argv[currentArg],"-title") == 0) {
	      simulationInfo.setTitle(argv[currentArg+1]);	
	      currentArg += 2;
	      
	    } else if (strcmp(argv[currentArg],"-description") == 0) {
	      simulationInfo.setDescription(argv[currentArg+1]);	
	      currentArg += 2;      
	    } else
	      more = false;
	  }
	} else 
	  currentArg++; 	
      }

      Tcl_SourceRCFile(interp);
      
      /*
       * Process commands from stdin until there's an end-of-file.  Note
       * that we need to fetch the standard channels again after every
       * eval, since they may have been changed.
       */
     
      if (simulationInfoOutputFilename != 0) {
	simulationInfo.start();
      }
     
      commandPtr = Tcl_NewObj();
      Tcl_IncrRefCount(commandPtr);
      
      inChannel = Tcl_GetStdChannel(TCL_STDIN);
      outChannel = Tcl_GetStdChannel(TCL_STDOUT);
      gotPartial = 0;
      while (1) {
	if (tty) {
	  Tcl_Obj *promptCmdPtr;
	  
	  char one[12] = "tcl_prompt1";
	  char two[12] = "tcl_prompt2";
	  promptCmdPtr = Tcl_GetVar2Ex(interp,
				       (gotPartial ? one : two),
				       NULL, TCL_GLOBAL_ONLY);
	  if (promptCmdPtr == NULL) {
	  defaultPrompt:
	    if (!gotPartial && outChannel) {
	      Tcl_WriteChars(outChannel, "OpenSees > ", 11);
	    }
	  } else {
	    
	    code = Tcl_EvalObjEx(interp, promptCmdPtr, 0);
	    
	    inChannel = Tcl_GetStdChannel(TCL_STDIN);
	    outChannel = Tcl_GetStdChannel(TCL_STDOUT);
	    errChannel = Tcl_GetStdChannel(TCL_STDERR);
	    if (code != TCL_OK) {
	      if (errChannel) {
		Tcl_WriteObj(errChannel, Tcl_GetObjResult(interp));
		Tcl_WriteChars(errChannel, "\n", 1);
	      }
	      Tcl_AddErrorInfo(interp,
			       "\n    (script that generates prompt)");
	      goto defaultPrompt;
	    }
	  }
	  if (outChannel) {
	    Tcl_Flush(outChannel);
	    }
	}
	if (!inChannel) {
	  goto done;
	}
	length = Tcl_GetsObj(inChannel, commandPtr);
	if (length < 0) {
	  goto done;
	}
	if ((length == 0) && Tcl_Eof(inChannel) && (!gotPartial)) {
	  goto done;
	}
	
	/*
	 * Add the newline removed by Tcl_GetsObj back to the string.
	 */
	
	Tcl_AppendToObj(commandPtr, "\n", 1);
	if (!TclObjCommandComplete(commandPtr)) {
	  gotPartial = 1;
	  continue;
	}
	
	gotPartial = 0;
	code = Tcl_RecordAndEvalObj(interp, commandPtr, 0);
	inChannel = Tcl_GetStdChannel(TCL_STDIN);
	outChannel = Tcl_GetStdChannel(TCL_STDOUT);
	errChannel = Tcl_GetStdChannel(TCL_STDERR);
	Tcl_DecrRefCount(commandPtr);
	commandPtr = Tcl_NewObj();
	Tcl_IncrRefCount(commandPtr);
	if (code != TCL_OK) {
	  if (errChannel) {
	    Tcl_WriteObj(errChannel, Tcl_GetObjResult(interp));
	    Tcl_WriteChars(errChannel, "\n", 1);
	  }
	} else if (tty) {
	  resultPtr = Tcl_GetObjResult(interp);
	  Tcl_GetStringFromObj(resultPtr, &length);
	  if ((length > 0) && outChannel) {
	    Tcl_WriteObj(outChannel, resultPtr);
	    Tcl_WriteChars(outChannel, "\n", 1);
	  }
	}
#ifdef TCL_MEM_DEBUG
	if (tclMemDumpFileName != NULL) {
	  Tcl_DecrRefCount(commandPtr);
	  Tcl_DeleteInterp(interp);
	  Tcl_Exit(0);
	}
#endif
      }
    }

 done:
    
    if (commandPtr != NULL) {
      Tcl_DecrRefCount(commandPtr);
    }
    
    
#ifdef _PARALLEL_PROCESSING
    return;
#endif
    
#ifdef _PARALLEL_INTERPRETERS
    return;
#endif
    
    /*
     * Rather than calling exit, invoke the "exit" command so that
     * users can replace "exit" with some other command to do additional
     * cleanup on exit.  The Tcl_Eval call should never return.
     */
    
    sprintf(buffer, "exit %d", exitCode);
    Tcl_Eval(interp, buffer);
    
    return;
}


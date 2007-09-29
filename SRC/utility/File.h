/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.1 $
// $Date: 2007-09-29 01:59:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/utility/File.h,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 09/07
//
// Description: This file contains the class definition for File used in SimulationINformation.
//
// What: "@(#) SimulationInformation.h, revA"

#ifndef File_h
#define File_h

#include <map>
#include <string>
using namespace std;

class FileIter;

class File
{
 public:
  File(const char *name, const char *description, bool isDir);
  ~File();

  int clear(void);
  int addFile(File *);
  const char *getName(void);
  const char *getDescription(void);
  bool isDir(void);
  FileIter getFiles(void);

  friend class FileIter;

 private:
  bool isDirectory;
  string name;
  string description;
  map<string, File *> dirFiles; 
};

#endif

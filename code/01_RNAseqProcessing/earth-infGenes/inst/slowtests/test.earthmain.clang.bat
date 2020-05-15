@rem test.earthmain.clang.bat: test the standalone earth.c with main()
@rem
@rem Stephen Milborrow Dec 2014 Shrewsbury

@cp "d:/bin/R320dll/i386/R.dll" .
                                @if %errorlevel% neq 0 goto error
@cp "d:/bin/R320dll/i386/Rblas.dll" .
                                @if %errorlevel% neq 0 goto error
@cp "d:/bin/R320dll/i386/Riconv.dll" .
                                @if %errorlevel% neq 0 goto error
@cp "d:/bin/R320dll/i386/Rgraphapp.dll" .
                                @if %errorlevel% neq 0 goto error
@cp "d:/bin/R320dll/i386/Rzlib.dll" .
                                @if %errorlevel% neq 0 goto error
@rem you may have to create Rdll.lib and Rblas.lib beforehand
@cp "../../.#/Rdll.lib" .
                                @if %errorlevel% neq 0 goto error
@cp "../../.#/Rblas.lib" .
                                @if %errorlevel% neq 0 goto error
@rem get iconv.dll from /a/r/ra/src/gnuwin32/unicode
@cp "../../.#/Rdll.lib" .
                                @if %errorlevel% neq 0 goto error

@rem modify the path to include clang, if needed
@set | egrep -i "^PATH=.*LLVM" >NUL && goto donesetpath
@echo Modifying path for clang
@set path=C:\Program Files (x86)\LLVM\bin;%PATH%
:donesetpath

@rem flags same as gcc, except 2nd and 3rd lines (starting with -Weverything) are new here
clang -DSTANDALONE -DMAIN -Wall -pedantic -Wextra -Weverything -O3 -std=gnu99^
 -Wno-missing-noreturn -Wno-float-equal -Wno-format-nonliteral^
 -Wno-padded -Wno-shadow -Wno-sign-conversion -Wno-undef^
 -I"/a/r/ra/include" -I../../inst/slowtests ../../src/earth.c^
 Rdll.lib Rblas.lib -o earthmain-clang.exe
                                @if %errorlevel% neq 0 goto error
@earthmain-clang.exe > test.earthmain-clang.out
                                @if %errorlevel% neq 0 goto error

@rem we use -w on mks.diff so it treats \r\n the same as \n
mks.diff -w test.earthmain-clang.out test.earthmain.out.save
                                @if %errorlevel% neq 0 goto error

@rm -f R.dll Rblas.dll Riconv.dll Riconv.dll Rgraphapp.dll Rzlib.dll Rdll.lib Rblas.lib earthmain-clang.* test.earthmain-clang.* *.o
@exit /B 0

:error
@exit /B %errorlevel%

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    blastlogic.cpp \
    entrezdialog.cpp \
    fileselectiondialog.cpp \
    main.cpp \
    mainwindow.cpp \
    optionsdata.cpp \
    optionsdialog.cpp \
    sequencedialog.cpp

HEADERS += \
    blastlogic.h \
    entrezdialog.h \
    fileselectiondialog.h \
    mainwindow.h \
    optionsdata.h \
    optionsdialog.h \
    sequencedialog.h

FORMS += \
    entrezdialog.ui \
    fileselectiondialog.ui \
    mainwindow.ui \
    optionsdialog.ui \
    sequencedialog.ui

QMAKE_CXXFLAGS += -I $$PWD/include
QMAKE_CXXFLAGS += -I $$PWD/include/ncbi-tools++
QMAKE_CXXFLAGS += -std=c++17
QMAKE_CXXFLAGS += -Wall
QMAKE_CXXFLAGS += -Wpedantic
QMAKE_CXXFLAGS += -save-temps=obj
QMAKE_CXXFLAGS += -O2
QMAKE_CXXFLAGS += -g3
QMAKE_CXXFLAGS += -pthread

QMAKE_LFLAGS += -L $$PWD/lib
QMAKE_LFLAGS += -Wl,-Map=FirstTest.map

LIBS += -lblastinput-static
LIBS += -lncbi_xloader_blastdb_rmt-static
LIBS += -lncbi_xloader_blastdb-static
LIBS += -lxblastformat-static
LIBS += -lalign_format-static
LIBS += -ltaxon1-static
LIBS += -lblastdb_format-static
LIBS += -lgene_info-static
LIBS += -lxformat-static
LIBS += -lxcleanup-static
LIBS += -lgbseq-static
LIBS += -lxobjedit-static
LIBS += -lefetch-static
LIBS += -leutils-static
LIBS += -legquery-static
LIBS += -lelink-static
LIBS += -lepost-static
LIBS += -lesearch-static
LIBS += -lespell-static
LIBS += -lesummary-static
LIBS += -leinfo-static
LIBS += -luilist-static
LIBS += -lehistory-static
LIBS += -lxobjread-static
LIBS += -lvariation-static
LIBS += -lsubmit-static
LIBS += -lxlogging-static
LIBS += -ltaxon3-static
LIBS += -lmlacli-static
LIBS += -lmla-static
LIBS += -lmedlars-static
LIBS += -lpubmed-static
LIBS += -lmedline-static
LIBS += -lbiblio-static
LIBS += -lgeneral-static
LIBS += -lvalid-static
LIBS += -lxalnmgr-static
LIBS += -lblastxml-static
LIBS += -lblastxml2-static
LIBS += -lxcgi-static
LIBS += -lxhtml-static
LIBS += -lproteinkmer-static
LIBS += -lxblast-static
LIBS += -lutrtprof-static
LIBS += -lxalgoblastdbindex-static
LIBS += -lcomposition_adjustment-static
LIBS += -lxalgodustmask-static
LIBS += -lxalgowinmask-static
LIBS += -lseqmasks_io-static
LIBS += -lseqdb-static
LIBS += -lblast_services-static
LIBS += -lxalnmgr-static
LIBS += -lxobjutil-static
LIBS += -lxobjread-static
LIBS += -lvariation-static
LIBS += -lsubmit-static
LIBS += -lxlogging-static
LIBS += -lxnetblastcli-static
LIBS += -lxnetblast-static
LIBS += -lblastdb-static
LIBS += -lscoremat-static
LIBS += -ltables-static
LIBS += -llmdb-static
LIBS += -lncbi_xloader_genbank-static
LIBS += -lncbi_xreader_id1-static
LIBS += -lncbi_xreader_id2-static
LIBS += -lncbi_xreader_cache-static
LIBS += -lncbi_xreader_pubseqos-static
LIBS += -lncbi_xreader_pubseqos2-static
LIBS += -ldbapi_driver-static
LIBS += -lncbi_xreader-static
LIBS += -lxconnect-static
LIBS += -lid1-static
LIBS += -lid2-static
LIBS += -lxobjmgr-static
LIBS += -lgenome_collection-static
LIBS += -lseqedit-static
LIBS += -lseqsplit-static
LIBS += -lsubmit-static
LIBS += -lseqset-static
LIBS += -lseq-static
LIBS += -lseqcode-static
LIBS += -lsequtil-static
LIBS += -lpub-static
LIBS += -lmedline-static
LIBS += -lbiblio-static
LIBS += -lgeneral-static
LIBS += -lxser-static
LIBS += -lxutil-static
LIBS += -lxncbi-static
LIBS += -lxcompress-static
LIBS += -lbz2-static
LIBS += -lz
LIBS += -lresolv

LIBS += -lxobjutil
LIBS += -lncbi_xloader_genbank
LIBS += -lncbi_xreader_id1
LIBS += -lncbi_xreader_id2
LIBS += -lncbi_xreader_cache
LIBS += -lncbi_xreader_pubseqos
LIBS += -lncbi_xreader_pubseqos2
LIBS += -ldbapi_driver
LIBS += -lncbi_xreader
LIBS += -lxconnect
LIBS += -lid1
LIBS += -lid2
LIBS += -lxobjmgr
LIBS += -lgenome_collection
LIBS += -lseqedit
LIBS += -lseqsplit
LIBS += -lsubmit
LIBS += -lseqset
LIBS += -lseq
LIBS += -lseqcode
LIBS += -lsequtil
LIBS += -lpub
LIBS += -lmedline
LIBS += -lbiblio
LIBS += -lgeneral
LIBS += -lxser
LIBS += -lxutil
LIBS += -lxncbi
LIBS += -lxcompress
LIBS += -lbz2

LIBS += -lpthread
LIBS += -ldl
LIBS += -lm

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

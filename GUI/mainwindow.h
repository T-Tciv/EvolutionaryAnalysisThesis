#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QListWidget>
#include "blastlogic.h"
#include "fileselectiondialog.h"
#include "optionsdialog.h"
#include "entrezdialog.h"
#include "sequencedialog.h"
#include <iostream>
#include <fstream>
#include <streambuf>
#include <string>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_fileSelectionButton_clicked();

    void on_listWidget_itemDoubleClicked(QListWidgetItem *item);

    void getFile(QString fileName);

    void getEntrez(QString entrez);

    void getOptions(OptionsData optionsSet);

    void on_pushButton_clicked();

    void on_inputRadioButton_toggled(bool checked);

    void on_fileRadioButton_toggled(bool checked);

    void on_filterButton_clicked();

    void on_searchButton_clicked();

private:
    Ui::MainWindow *ui;
    BlastLogic *blastLogic;
    FileSelectionDialog *fileDialog;
    OptionsDialog *optionsDialog;
    EntrezDialog *entrezDialog;
    std::vector<std::string> blastVector;
};
#endif // MAINWINDOW_H

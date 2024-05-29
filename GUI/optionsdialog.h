#ifndef OPTIONSDIALOG_H
#define OPTIONSDIALOG_H

#include <QDialog>
#include <QMenu>
#include <QList>
#include <QAction>
#include "optionsdata.h"

namespace Ui {
class OptionsDialog;
}

class OptionsDialog : public QDialog
{
    Q_OBJECT

public:
    explicit OptionsDialog(QWidget *parent = nullptr);
    // TODO: ДОБАВИТЬ КОНСТРУКТОР ДЛЯ УСТАНОВКИ ДОПУСТИМЫХ ЗНАЧЕНИЙ ПРИ СОЗДАНИИ ОКНА ИЗВНЕ
    ~OptionsDialog();

private slots:
    void on_setOptionsButton_clicked();

signals:
    void setOptions(OptionsData optionsSet);

private:
    Ui::OptionsDialog *ui;

    OptionsData *optionsData;

    QMenu *programMenu;
    QList<QAction> *programList;

    QMenu *databaseMenu;
    QList<QAction> *databaseList;
};

#endif // OPTIONSDIALOG_H

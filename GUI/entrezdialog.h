#ifndef ENTREZDIALOG_H
#define ENTREZDIALOG_H

#include <QDialog>
#include <iostream>
#include <string>

namespace Ui {
class EntrezDialog;
}

class EntrezDialog : public QDialog
{
    Q_OBJECT

public:
    explicit EntrezDialog(QWidget *parent = nullptr);
    ~EntrezDialog();

private slots:
    void on_setEntrezButton_clicked();

signals:
    void setEntrez(QString entrezQuery);

private:
    Ui::EntrezDialog *ui;
};

#endif // ENTREZDIALOG_H

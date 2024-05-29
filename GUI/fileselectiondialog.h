#ifndef FILESELECTIONDIALOG_H
#define FILESELECTIONDIALOG_H

#include <QDialog>
#include <QFileSystemModel>
#include <QAbstractItemModel>
#include <iostream>
#include <string>
#include <memory>

namespace Ui {
class FileSelectionDialog;
}

class FileSelectionDialog : public QDialog
{
    Q_OBJECT

public:
    explicit FileSelectionDialog(QWidget *parent = nullptr);
    ~FileSelectionDialog();
    std::string getFileName();

private slots:
    void on_listView_doubleClicked(const QModelIndex &index);

    void on_fileSelectionButton_clicked();

signals:
    void setFile(QString fileName);

private:
    Ui::FileSelectionDialog *ui;
    std::unique_ptr<QFileSystemModel> model;
};

#endif // FILESELECTIONDIALOG_H

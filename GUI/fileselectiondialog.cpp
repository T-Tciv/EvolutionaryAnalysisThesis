#include "fileselectiondialog.h"
#include "ui_fileselectiondialog.h"

FileSelectionDialog::FileSelectionDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::FileSelectionDialog)
{
    ui->setupUi(this);
    model = std::make_unique<QFileSystemModel>(new QFileSystemModel(this));
    model->setNameFilters(QStringList() << "*.fst");

    //TODO: возможно, на Windows тут ещё QDir:Drives нужен, но на линуксе он игнорируется, если верить документации
    model->setFilter(QDir::AllDirs | QDir::Files | QDir::NoDot);

    model->setRootPath(QDir::rootPath());
    ui->listView->setModel(model.get());
    ui->listView->setRootIndex(model->index(QDir::homePath()));
}

FileSelectionDialog::~FileSelectionDialog()
{
    delete ui;
}

std::string FileSelectionDialog::getFileName() {
    return "File name sample";
}

void FileSelectionDialog::on_listView_doubleClicked(const QModelIndex &index)
{
    // Приводим объект-источник сигнала к типу QListView
    QListView *  listView = (QListView*) sender();
    // С помощью модели получаем данные о файле
    QFileInfo fileInfo = model->fileInfo(index);

    // Теперь используем данные о файле, чтобы обработать нажатие
    QString fileName = fileInfo.fileName();
    bool isDir = fileInfo.isDir();

    if ("." == fileName) {
        listView->setRootIndex(model->index(""));
    } else if (".." == fileName) {
        QDir dir = fileInfo.dir();
        dir.cdUp();
        QModelIndex rootIndex = model->index(dir.absolutePath());
        listView->setRootIndex(rootIndex);
    } else if (isDir) {
        listView->setRootIndex(index);
    }
}

void FileSelectionDialog::on_fileSelectionButton_clicked()
{
    QFileInfo fileInfo = model->fileInfo(ui->listView->currentIndex());
    QString fileName = fileInfo.absoluteFilePath();

    // Выбирать разрешаем только файлы, папки не трогаем
    if (fileInfo.isFile()) {
        emit setFile(fileName);
        close();
    }

}


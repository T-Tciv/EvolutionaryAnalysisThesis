#include "entrezdialog.h"
#include "ui_entrezdialog.h"

EntrezDialog::EntrezDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::EntrezDialog)
{
    ui->setupUi(this);
}

EntrezDialog::~EntrezDialog()
{
    delete ui;
}

void EntrezDialog::on_setEntrezButton_clicked()
{
    QString query = "";
    std::vector<QString> entrezVector;

    QString organismNameText = ui->organismField->text();
    QString organismName = organismNameText + "[Organism]";

    QString propStringText = ui->typeField->text();
    QString props = propStringText + "[Properties]";

    QString keyWordsText = ui->keyWordField->text();
    QString keyWords = keyWordsText + "[Keyword]";

    if (!organismNameText.isEmpty())
        entrezVector.push_back(organismName);

    if (!propStringText.isEmpty())
        entrezVector.push_back(props);

    if (!keyWordsText.isEmpty())
        entrezVector.push_back(keyWords);

    if (entrezVector.size() != 0) {
        query += entrezVector.back();
        entrezVector.pop_back();
    }

    for (QString entrez : entrezVector) {
        query += "AND " + entrez;
    }

    emit setEntrez(query);
    close();
}


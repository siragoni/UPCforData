{
 gROOT->Reset();
 AliCDBManager *man = AliCDBManager::Instance();
 man->Init();
 //man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
 man->SetDefaultStorage("alien://folder=/alice/data/2009/OCDB/");
 //1
 man->SetRun(104065);
 AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/CTP/Config");
 AliTriggerConfiguration* rc = dynamic_cast<AliTriggerConfiguration*>(entry->GetObject());
 rc->Print();
 //2
 // man->SetRun(104066);
 // AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/CTP/Config");
 // AliTriggerConfiguration* rc = dynamic_cast<AliTriggerConfiguration*>entry->GetObject();
 // rc->Print();
 // //3
 // man->SetRun(104067);
 // AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/CTP/Config");
 // AliTriggerConfiguration* rc = dynamic_cast<AliTriggerConfiguration*>entry->GetObject();
 // rc->Print();
 // //4
 // man->SetRun(104068);
 // AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/CTP/Config");
 // AliTriggerConfiguration* rc = dynamic_cast<AliTriggerConfiguration*>entry->GetObject();
 // rc->Print();
 // //5
 // man->SetRun(104070);
 // AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/CTP/Config");
 // AliTriggerConfiguration* rc = dynamic_cast<AliTriggerConfiguration*>entry->GetObject();
 // rc->Print();
 // //6
 // man->SetRun(104071);
 // AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/CTP/Config");
 // AliTriggerConfiguration* rc = dynamic_cast<AliTriggerConfiguration*>entry->GetObject();
 // rc->Print();
 // //7
 // man->SetRun(104073);
 // AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/CTP/Config");
 // AliTriggerConfiguration* rc = dynamic_cast<AliTriggerConfiguration*>entry->GetObject();
 // rc->Print();
 // //8
 // man->SetRun(104080);
 // AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/CTP/Config");
 // AliTriggerConfiguration* rc = dynamic_cast<AliTriggerConfiguration*>entry->GetObject();
 // rc->Print();
 // //9
 // man->SetRun(104083);
 // AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/CTP/Config");
 // AliTriggerConfiguration* rc = dynamic_cast<AliTriggerConfiguration*>entry->GetObject();
 // rc->Print();
 // //10
 // man->SetRun(104152);
 // AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/CTP/Config");
 // AliTriggerConfiguration* rc = dynamic_cast<AliTriggerConfiguration*>entry->GetObject();
 // rc->Print();
 // //11
 // man->SetRun(104155);
 // AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/CTP/Config");
 // AliTriggerConfiguration* rc = dynamic_cast<AliTriggerConfiguration*>entry->GetObject();
 // rc->Print();
 // //12
 // man->SetRun(104157);
 // AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/CTP/Config");
 // AliTriggerConfiguration* rc = dynamic_cast<AliTriggerConfiguration*>entry->GetObject();
 // rc->Print();
 // //13
 // man->SetRun(104158);
 // AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/CTP/Config");
 // AliTriggerConfiguration* rc = dynamic_cast<AliTriggerConfiguration*>entry->GetObject();
 // rc->Print();
 // //14
 // man->SetRun(104159);
 // AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/CTP/Config");
 // AliTriggerConfiguration* rc = dynamic_cast<AliTriggerConfiguration*>entry->GetObject();
 // rc->Print();

}

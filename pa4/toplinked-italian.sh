#!/bin/bash

sbt package

/usr/local/spark-2.1.0-bin-hadoop2.4/bin/spark-submit \
        --class TopLinked\
        --master local\
        --deploy-mode client\
        --packages com.databricks:spark-xml_2.11:0.4.1 \
        `pwd`/target/scala-2.11/cs149-spark-assignment_2.11-1.0.jar \
                `cat /usr/local/etc/master`\
                s3n://cs149-spark-pa4/itwiki-20170120-pages-articles-multistream.xml


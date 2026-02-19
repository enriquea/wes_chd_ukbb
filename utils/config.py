"""
Environment-specific path configuration.

Configure these paths for your environment by setting environment variables
before running any pipeline script:

    export WES_NFS_DIR="file:///mnt/mydata"
    export WES_HDFS_DIR="hdfs://my-cluster:9820"

If the variables are not set, the original production defaults are used so
that existing behaviour is preserved.
"""

import os

NFS_DIR = os.environ.get("WES_NFS_DIR", "file:///home/ubuntu/data")
HDFS_DIR = os.environ.get("WES_HDFS_DIR", "hdfs://spark-master:9820")

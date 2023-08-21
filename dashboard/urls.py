from django.urls import path
from .views import *

app_name = "dashboard"

urlpatterns = [
    path('users-dashboard/<user>/', DashBoardView.as_view(), name="anasayfa"),
    path('user-posts/<user>/<int:pk>/', MyPostList.as_view(), name="user_posts"),
    path('profile/<user>/', ProfileListView.as_view(), name="profile_dash"),
    path('social_media/<user>/', SocialMediaListView.as_view(), name="social_dash"),
    path('messages/<user>/', MessagesListView.as_view(), name="messages_list"),
    path('profile-update/<username>/<pk>/', profile_update, name="profile_update"),
    path('social_media/<user>/add/', SocialMediaCreateView.as_view(), name="social_create"),
    path('social_update/<user>/<pk>/', SocialMediaUpdateView.as_view(), name="social_update"),
    path('messages_detail/<sender>/<pk>/', MessagesDetailView.as_view(), name="message_detail"),
    path('read/<pk>/', mark_read_message, name="mark_read_message"),
    path('mark_as_read_all/', mark_as_read_all, name="mark_as_read_all"),
    path('delete_message/<pk>/', message_delete, name="message_delete"),
    path('report_message/<pk>/', report_message, name="report_message"),
    path('user-reply-message/<int:pk>/<username>/<int:user_pk>/', user_reply_message, name="user_reply_message"),
    path('user-sent-message/<int:pk>/<username>/', user_sent_message, name="user_sent_message"),
]
